#include <gtest/gtest.h>

#include "stored.h"
#include "sciplot/sciplot.hpp"
#include "gnuplot-iostream.h"

TEST(test5g, test5g)
{
    const char *filePathEnv = std::getenv("FILE_PATH");
    const char *hdfNameEnv = std::getenv("HDF_NAME");
    std::string homeDir = getenv("HOME");
    std::filesystem::path backOneDir = std::filesystem::current_path();
    std::string filePath = filePathEnv ? homeDir + filePathEnv : "../IQbin/IQbin/UE1-0503-1857.dat";
    std::string hdfName = hdfNameEnv ? hdfNameEnv : "default";
    stftTest testSignal;
    std::vector<double> window;
    std::vector<std::complex<double>> iqBuff;
    double fs;
    double fDelta;
    double tDelta;
    std::streampos begin, end;
    bool testing = true;

    std::vector<double> inphaseOrig;
    std::vector<double> quadOrig;
    std::vector<double> i2;
    std::vector<double> q2;
    std::vector<double> signalRecreated;
    std::vector<double> signalRecreated2;
    // std::ofstream outfile("UE1-0503-1857-test3.csv", std::ios::out | std::ios::trunc);
    // std::ofstream testfile("testfile.csv", std::ios::out | std::ios::trunc);
    std::ifstream file(filePath, std::ios::binary | std::ios::ate | std::ios::in);

    if (!file.is_open())
    {
        std::cout << "Error opening file" << std::endl;
        std::cout << "File path: " << filePath << std::endl;
    }

    std::streampos size = file.tellg();
    std::cout << "Size of file: " << size << "\n";
    if (size == -1)
    {
        std::cout << "Error getting file size" << std::endl;
    }
    int totalSample = size / (2 * sizeof(int16_t));
    unsigned int chunkSize = 768000 * (2 * sizeof(int16_t));
    const int seeker = 1;
    int count = 0;
    file.seekg(8, std::ios::beg);

    std::vector<unsigned char> buffer(chunkSize);
    file.read((char *)buffer.data(), chunkSize);

    for (int i = 0; i < buffer.size(); i = i + 2)
    {
        iqBuff.emplace_back(std::complex<float>(buffer[i], buffer[i + 1]));
    }

    // file.seekg(chunkSize, std::ios::cur);
    std::cout << "Samples : " << iqBuff.size() << "\n";

    sfft(window, testSignal, iqBuff);

    // saving to hdf5 file

    fs = 40e6;
    fDelta = floor(fs / testSignal.windowSize);
    tDelta = testSignal.windowSize / fs;

    std::cout << "Ifft Output Size: " << testSignal.ifft.size() << std::endl;
    for (auto &z : testSignal.ifft)
    {
        double sum;
        sum = std::abs(z);
        i2.push_back(std::real(z) / sum);
        q2.push_back(std::imag(z) / sum);
    }
    int index = 0;
  

    for (auto &n : iqBuff)
    {
        double sum;
        sum = std::abs(n);
        inphaseOrig.push_back(std::real(n) / sum);
        quadOrig.push_back(std::imag(n) / sum);
    }

  

    std::complex<double> mseCom = 0;
    double mse;
    double min, max;
    std::vector<std::complex<double>> diff;

    for (int i = 0; i < iqBuff.size(); i++)
    {
        diff.push_back(iqBuff[i] - testSignal.ifft[i]);
        mseCom += diff[i];

    }
    mse = std::abs(mseCom) / iqBuff.size();
  
  

    std::vector<double> xDecimated;
    std::vector<double> xDecimated2;
    std::vector<double> yDecimated;
    std::vector<double> yDecimated2;
    std::vector<double> difference;
    int N = testSignal.ifft.size();
    sciplot::Vec x = sciplot::linspace(0, N / fs, iqBuff.size());
    int decimationFactor = 100;
    for (size_t i = 0; i < iqBuff.size(); i += decimationFactor)
    {
        xDecimated.push_back(inphaseOrig[i]);
        yDecimated.push_back(quadOrig[i]);
        xDecimated2.push_back(i2[i]);
        yDecimated2.push_back(q2[i]);
        signalRecreated.push_back(inphaseOrig[i] * std::cos(2 * M_PI * i) + quadOrig[i] * std::sin(2 * M_PI * i));
        signalRecreated2.push_back(i2[i] * std::cos(2 * M_PI * i) + q2[i] * std::sin(2 * M_PI * i));
    }
    for(int i = 0; i < N; i+=decimationFactor)
    {
        difference.push_back(std::abs(diff[i]));

    }
    std::vector<double> normDiff;
    min = *std::min_element(difference.begin(), difference.end());
    max = *std::max_element(difference.begin(), difference.end());
    for (auto &val : difference)
    {
        normDiff.push_back((val-min) / (max-min));
    }
    std::cout << "MSE: " << mse << std::endl;
    sciplot::Plot2D plot;
    plot.size(1400, 800);
    plot.xlabel("TIME");
    plot.ylabel("Percentage Difference");

    plot.autoclean(false);
    plot.drawCurve(x, normDiff).label("difference");
    

    sciplot::Figure fig = {{plot}};
    sciplot::Canvas canvas = {{fig}};
    canvas.size(1080, 600);
    canvas.save("test-qam.pdf");

    std::cout << "Frequency Resolution: " << fDelta / 1000 << " kHz \n";
    std::cout << "Window Duration: " << tDelta * 1000 << " ms \n";
    std::cout << "Total Time: " << tDelta + (testSignal.windowCount - 1) * (testSignal.hopSize / fs) * 1000 << " ms \n";
    
    EXPECT_NEAR(mse, 0.0, 0.01);

    file.close();

    std::vector<std::vector<double>> abs;

    for (auto &row : testSignal.stftAbs)
    {
        double min = *std::min_element(row.begin(), row.end());
        double max = *std::max_element(row.begin(), row.end());
        std::vector<double> absRow;
        for (auto &value : row)
        {
            if ((value - min) == 0 || (max - min) == 0)
            {
                absRow.push_back(-100);
            }
            else
            {
                absRow.push_back(10 * (log10(((value - min))) - log10((max - min))));
            }
        }
        abs.push_back(absRow);
    }
    std::ofstream outFile(hdfName + ".txt", std::ios::trunc);
    if (!outFile)
    {
        std::cerr << "Error opening file for writing: " << hdfName << std::endl;
        return;
    }

    for (const auto &row : abs)
    {
        for (const auto &value : row)
        {
            outFile << value << " ";
        }
        outFile << "\n";
    }

    outFile.close();
    std::cout << "Data exported to " << hdfName << std::endl;
    Gnuplot gp;

    // Set Gnuplot parameters
    gp << "set terminal pngcairo size 1080,600\n";
    gp << "set output 'spectrogram2.png'\n";
    gp << "set pm3d map\n";
    gp << "set palette defined (0 'blue', 1 'cyan', 2 'green', 3 'yellow', 4 'red')\n";
    gp << "set xlabel 'Time'\n";
    gp << "set ylabel 'Frequency'\n";
    gp << "set cblabel 'Magnitude'\n";
    gp << "splot 'testWave.txt'  matrix using 2:1:3 with image\n";
}

int main(int argc, char **argv)
{

    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
