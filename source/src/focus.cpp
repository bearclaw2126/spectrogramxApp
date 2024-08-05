#include "stored.h"
#include <filesystem>
#include "sciplot/sciplot.hpp"
int main(int argc, char *argv[])
{
   
    const char* filePathEnv = std::getenv("FILE_PATH");
    const char* hdfNameEnv = std::getenv("HDF_NAME");
    std::string homeDir = getenv("HOME");
    std::filesystem::path backOneDir = std::filesystem::current_path();
    std::string filePath = filePathEnv ? homeDir +filePathEnv : "../IQbin/IQbin/UE1-0503-1857.dat";
    std::string hdfName = hdfNameEnv ? hdfNameEnv : "default";
    stftTest sig;
    std::vector<double> window;
    std::vector<std::complex<double>> iqBuff;
   
    double fs;
    double fDelta;
    double tDelta;
    std::streampos begin, end;

    
    // std::ofstream outfile("UE1-0503-1857-test3.csv", std::ios::out | std::ios::trunc);
    // std::ofstream testfile("testfile.csv", std::ios::out | std::ios::trunc);
    std::ifstream file(filePath, std::ios::binary  | std::ios::ate | std::ios::in);
    
    if (!file.is_open())
    {
        std::cout << "Error opening file" << std::endl;
        std::cout << "File path: " << filePath << std::endl;
        return -1;
    }

    std::streampos size = file.tellg();
    std::cout << "Size of file: " << size << "\n";
   if (size == -1)
   {
         std::cout << "Error getting file size" << std::endl;
   }
    int totalSample = size / (2 * sizeof(int16_t));
    unsigned int chunkSize = 4096 * (2 * sizeof(int16_t));
    const int seeker = 100;
    int count = 0;
    file.seekg(8, std::ios::beg); 

    for (int i = 0; i < seeker; ++i)
    {

        std::vector<unsigned char> buffer(chunkSize);
        file.read((char *)buffer.data(), chunkSize);

        for (int i = 0; i < buffer.size(); i = i + 2)
        {
            iqBuff.emplace_back(std::complex<float>(buffer[i], buffer[i + 1]));
        }
        

        file.seekg(chunkSize, std::ios::cur);
        std::cout << "Samples : " << iqBuff.size() << "\n";
    
    
        sfft(window,sig, iqBuff);
       
    
        count++;


        // saving to hdf5 file
       IQfile(iqBuff, hdfName);
      
        //IQfile(sig, hdfName);
        fs = 40e6;
        fDelta = floor(40e6 / sig.stftAbs[0].size());
        tDelta = sig.stftAbs[0].size() / fs;
        // tDelta  = 1 / fDelta;
        
        sig.stftCIR.clear();
        sig.output.clear();
        iqBuff.clear();
      
    
    }
    

    std::cout << "Frequency Resolution: " << fDelta / 1000 << " kHz \n";
    std::cout << "Window Duration: " << tDelta * 1000 << " ms \n";
    std::cout << "Total Time: " << tDelta + (sig.windowCount - 1) * (sig.hopSize / fs) * 1000 << " ms \n";

  file.close();
    }

