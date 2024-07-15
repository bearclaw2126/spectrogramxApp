
#include "stored.h"
#include "H5Cpp.h"
void hamming(int w, std::vector<double> &window)
{
    for (int i = 0; i < w; ++i)
    {
        window.push_back(0.54 - 0.46 * cos(2 * M_PI * i / (w - 1)));
    }
}

stft stftCalc(std::vector<double> &window, std::vector<std::complex<double>> &iq)
{ // setting time counters to measure fft time
    using std::chrono::duration;
    using std::chrono::duration_cast;
    using std::chrono::high_resolution_clock;
    using std::chrono::milliseconds;
    auto t1 = high_resolution_clock::now();
    stft signalIn;
    int fs = 40e6;
    int scs = 30e3;
    int index = log2(ceil(fs / (scs * 12)));
    int signalLength = iq.size();
    signalIn.windowSize = windowSizeLUT[2];
    int windowSize = signalIn.windowSize;
    int hopSize = windowSize / 2;
    signalIn.hopSize = hopSize;
    int N = windowSize;
    int windowCount = floor((signalLength - windowSize) / hopSize);
    signalIn.windowCount = windowCount;
    fftw_plan p;

    std::cout << " WINDOW COUNT: " << windowCount << "\n";
    std::cout << " WINDOW SIZE: " << windowSize << "\n";
    std::cout << " HOP SIZE: " << hopSize << "\n";

    fftw_complex *in = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * windowSize * windowCount);
    fftw_complex *out = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * windowSize * windowCount);

    hamming(hopSize, window);

    int readingIndex = 0;

    for (int i = 0; i < windowCount; ++i)
    {
        int startIndex = i * hopSize;

        for (int j = 0; j < windowSize; ++j)
        {
            int signalIndex = startIndex + j;
            if (signalIndex < iq.size())
            {
                std::complex<double> windowSample = iq[signalIndex] * window[j];
                in[signalIndex][0] = std::real(windowSample);
                in[signalIndex][1] = std::imag(windowSample);
            }
            else
            {
                in[signalIndex][1] = 0;
                in[signalIndex][0] = 0;
            }
        }
    }

    std::cout << "Starting STFT\n";

    // int odist;
    // int ostride;
    // const int dim = 1;
    // const int n = windowSize;    // size of fft
    // int howMany = windowCount;   // number of windows
    // int idist = odist = windowSize; // determines the hop size, distance between the start of the transform to the next
    // int istride = ostride = 1;   // will take fft every set number of samples

    // p = fftw_plan_many_dft(dim, &n, howMany, in, NULL, istride, idist, out, NULL, ostride, odist, FFTW_FORWARD, FFTW_ESTIMATE);
    p = fftw_plan_dft_1d(windowSize, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

    signalIn.stftBuff.reserve(windowCount);
    for (int i = 0; i < windowCount * hopSize; i += hopSize)
    {
        fftw_execute(p);
        std::vector<std::complex<double>> temp(windowSize / 2 + 1);
        for (int j = 0; j < windowSize / 2 + 1; ++j)
        {
            int index = i + j;

            // temp.push_back(std::complex<double>(out[index][0], out[index][1]));
            temp[j] = std::complex<double>(out[j][0], out[j][1]);
            // Convert fftw_complex to std::complex
        }
        signalIn.stftBuff.push_back(temp);
    }

    std::vector<std::vector<std::complex<double>>>::iterator var;
    std::vector<std::complex<double>>::iterator win;

    auto t2 = high_resolution_clock::now();

    duration<double, std::milli> ms_double = t2 - t1;

    std::cout << "Time to execute FFT for 1 thread" << ": " << ms_double.count() << "ms \n";
    std::cout << "size of stft " << signalIn.stftBuff[0].size() << "\n";
    std::cout << "windows of stft " << signalIn.stftBuff.size() << "\n";
    // free memory
    fftw_destroy_plan(p);
    fftw_free(in);
    fftw_free(out);

    std::cout << "Calculating Abs\n";

    for (var = signalIn.stftBuff.begin(); var != signalIn.stftBuff.end(); ++var)
    {
        std::vector<double> temp;
        for (win = var->begin(); win != var->end(); ++win)
        {
            temp.push_back(20 * log10(std::abs(*win)) + .0000000001);
        }
        signalIn.stftAbs.push_back(temp);
    }
    std::cout << "size of stftAbs " << signalIn.stftAbs[0].size() << "\n";
    std::cout << "windows of stftAbs " << signalIn.stftAbs.size() << "\n";

    std::cout << "Calculated Abs\n";

    // return the newly created struct with the array

    return signalIn;
}

void cirValue(stft &s)
{

    std::cout << "Calculating CIR\n";
    // converting to log form to avoid having to deal with divsion first
    // CIR spectral power is repersented as 20*(log10(abs(signal[x+1] - log10(abs(signal[x]))
    int counter = 0;
    std::vector<std::vector<double>>::iterator var;
    std::vector<double>::iterator win;
    for (var = s.stftAbs.begin(); var != s.stftAbs.end(); ++var)
    {
        s.stftCIR.reserve(s.stftAbs.size());
        std::vector<double> temp;
        for (win = var->begin(); win != var->end(); ++win)
        {
            temp.push_back(10 * (*(win + 1) - *win));
        }
        s.stftCIR.push_back(temp);
    }
}

/**
 * @ingroup Functions
 * @brief This function writes the CIR values to an hdf5 file as a matrix array.
 * @param in
 */
void writeHdf5(stft &in, int round, std::string hdfName)
{

    // Define the directory and file name based on the round value
    std::filesystem::path backOneDir = std::filesystem::current_path().parent_path();

    std::string directory = (backOneDir / "database/").string();
    std::string fileName = "CIRUE1_" + std::to_string(round) + ".h5";
    std::string fullPath = directory + fileName;

    hid_t file_id, dataset_id, dataspace_id; // Identifiers for the first dataset
    herr_t status;
    hsize_t dims[2] = {in.stftCIR.size(), in.stftCIR[0].size()};
    size_t totalSize = in.stftCIR.size() * in.stftCIR[0].size();
    double *data = new double[totalSize];

    for (size_t i = 0; i < in.stftCIR.size(); ++i)
    {
        for (size_t j = 0; j < in.stftCIR[i].size(); ++j)
        {
            data[i * in.stftCIR[0].size() + j] = in.stftCIR[i][j];
        }
    }

    file_id = H5Fcreate(fullPath.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if (file_id < 0)
        goto cleanup_data;

    dataspace_id = H5Screate_simple(2, dims, NULL);
    if (dataspace_id < 0)
        goto cleanup_file;

    dataset_id = H5Dcreate2(file_id, "fft_output", H5T_IEEE_F64LE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (dataset_id < 0)
        goto cleanup_dataspace;

    status = H5Dwrite(dataset_id, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
    if (status < 0)
        goto cleanup_dataset;

cleanup_dataset:
    H5Dclose(dataset_id);
cleanup_dataspace:
    H5Sclose(dataspace_id);
cleanup_file:
    H5Fclose(file_id);
cleanup_data:
    delete[] data;
    ///////////////////////////////////////////////////////////////////////////////////////

    directory = (backOneDir / "./Fullchannel/").string(); // Specify the directory path
    fileName = "UE1STFT_" + std::to_string(round) + ".h5";
    fullPath = directory + fileName;
    data = new double[totalSize];
    // Fill stftAbsData from in.stftAbs

    for (size_t i = 0; i < in.stftCIR.size(); ++i)
    {
        for (size_t j = 0; j < in.stftCIR[i].size(); ++j)
        {
            data[i * in.stftAbs[0].size() + j] = in.stftAbs[i][j];
        }
    }

    file_id = H5Fcreate(fullPath.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if (file_id < 0)
        goto cleanup_data;

    dataspace_id = H5Screate_simple(2, dims, NULL);
    if (dataspace_id < 0)
        goto cleanup_file;

    dataset_id = H5Dcreate2(file_id, "fft_output", H5T_IEEE_F64LE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (dataset_id < 0)
        goto cleanup_dataspace;

    status = H5Dwrite(dataset_id, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
    if (status < 0)
        goto cleanup_dataset;
}

void signalHandler(double signal)
{
    if (signal == SIGSEGV)
    {
        std::cerr << "Error: Segmentation fault detected." << std::endl;
        // Perform any cleanup if necessary
        exit(1); // Exit the program with an error code
    }
}

void hdf5file(stft &signal, std::string fileName)
{
    std::filesystem::path backOneDir = std::filesystem::current_path().parent_path();
    std::string hdfName = fileName;
    bool test = false;
    std::string folder = test == true ? "/testOutput" : "database/";
    std::string directory = (backOneDir / folder).string();
    std::string fullPath = directory + hdfName + ".h5";
    hsize_t dims[2] = {signal.stftCIR.size(), signal.stftCIR[0].size()};
    std::cout << fullPath << "\n";
    H5::H5File file;
    if (!std::filesystem::exists(fullPath))
    {
        file = H5::H5File(fullPath, H5F_ACC_TRUNC);
        H5::Group group = file.createGroup("/CIR");
    }
    else
    {
        file = H5::H5File(fullPath, H5F_ACC_RDWR);
    }
    // Open the group and count the number of datasets
    H5::Group cir = file.openGroup("/CIR");
    int numDatasets = cir.getNumObjs();

    if (numDatasets <= 0)
    {
       std::string dataSetName = "CIR0";
        H5::DataSpace dataspace(2, dims);
        H5::DataSet dataset = cir.createDataSet(dataSetName, H5::PredType::NATIVE_DOUBLE, dataspace);
        dataset.write(signal.stftCIR.data(), H5::PredType::NATIVE_DOUBLE);
      
    }
    else if (numDatasets >= 1 && numDatasets <= 10)
    {
        std::string dataSetName = "CIR" + std::to_string(numDatasets);
        H5::DataSpace dataspace(2, dims);
        H5::DataSet dataset = cir.createDataSet(dataSetName, H5::PredType::NATIVE_DOUBLE, dataspace);
        dataset.write(signal.stftCIR.data(), H5::PredType::NATIVE_DOUBLE);
       
    }



    return;
}