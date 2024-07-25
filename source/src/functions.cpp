
#include "stored.h"
#include "H5Cpp.h"
void hamming(int w, std::vector<double> &window)
{
    for (int i = 0; i < w; ++i)
    {
        window.push_back(0.54 - 0.46 * cos(2 * M_PI * (i / (w - 1))));
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
    int signalLength = iq.size();
    signalIn.windowSize = windowSizeLUT[1];
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

    fftw_complex *in = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * windowSize);
    fftw_complex *out = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * windowSize);

    hamming(hopSize, window);

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

    int chunk = 0;
    int readIndex = 0;
    int bStop = 0;
    int windowIndex = 0;
    int chunkPosition = 0;

    // Allow for cases where signal is below and zero pad it
    std::vector<std::complex<double>> temp;
    temp.reserve(signalLength);
    for (int windowIndex = 0; windowIndex < windowCount; windowIndex++)
    {
        if (chunkPosition > signalLength - windowSize)
        {
            break;
        }
        else
        {
            for (int i = 0; i < windowSize; ++i)
            {
                readIndex = chunkPosition + i;

                if (readIndex < signalLength)
                {
                    in[i][0] = std::real(iq[readIndex]) * window[i];
                    in[i][1] = std::imag(iq[readIndex]) * window[i];
                }
                else
                {
                    in[i][0] = 0;
                    in[i][1] = 0;
                }
            }
            fftw_execute(p);

            for (int i = 0; i < windowSize; i++)
            {
                temp[chunkPosition + i] += std::complex<double>(out[i][0], out[i][1]);
            }

            chunkPosition += hopSize;
        }
    }
    signalIn.stftBuff.resize(windowCount);
    std::cout << "APPLYING TO STFT BUFF\n";
    // Reserve space for windowCount number of elements
    // signalIn.stftBuff[0].resize(hopSize);
    for (int i = 0; i < windowCount; i++)
    {
        for (int j = 0; j < windowSize; j++)
        {
            signalIn.stftBuff[i].push_back(temp[i * hopSize + j]);
        }
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
            temp.push_back(std::log10(std::abs(*win)) + .0000000001);
        }
        signalIn.stftAbs.push_back(temp);
    }
    std::cout << "size of stftAbs " << signalIn.stftAbs[0].size() << "\n";
    std::cout << "windows of stftAbs " << signalIn.stftAbs.size() << "\n";

    std::cout << "Calculated Abs\n";

    // return the newly created struct with the array

    return signalIn;
}

void cirValue(stftTest &s)
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
            temp.push_back(20 * (*(win + 1) - *win));
        }
        s.stftCIR.push_back(temp);
    }
    std::cout << "Calculated CIR\n"
              << s.stftCIR[0].size() << std::endl;
}

/**
 * @ingroup Functions
 * @brief This function writes the CIR values to an hdf5 file as a matrix array.
 * @param in
 */
void writeHdf5(stftTest &in, int round, std::string hdfName)
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

void hdf5file(stftTest &signal, std::string fileName)
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
}

void hdfSTFTtest(const stftTest signal, const std::string fileName)
{
    std::filesystem::path backOneDir = std::filesystem::current_path().parent_path();
    std::string hdfName = fileName;
    bool test = true;
    std::string folder = test == true ? "testOutput/" : "database/";
    std::string directory = (backOneDir / folder).string();
    std::string fullPath = directory + hdfName + ".h5";
    hsize_t dims[2] = {signal.stftAbs[0].size(), signal.stftAbs.size()};
    H5::H5File file;
    if (!std::filesystem::exists(fullPath))
    {
        file = H5::H5File(fullPath, H5F_ACC_TRUNC);
        H5::Group group = file.createGroup("/fftOutput");
    }
    else
    {
        file = H5::H5File(fullPath, H5F_ACC_RDWR);
    }
    // Open the group and count the number of datasets
    H5::Group cir = file.openGroup("/fftOutput");
    int numDatasets = cir.getNumObjs();

    if (numDatasets <= 0)
    {
        std::string dataSetName = "STFT0";
        H5::DataSpace dataspace(2, dims);
        H5::DataSet dataset = cir.createDataSet(dataSetName, H5::PredType::NATIVE_DOUBLE, dataspace);
        dataset.write(signal.stftAbs.data(), H5::PredType::NATIVE_DOUBLE);
    }
    else if (numDatasets >= 1 && numDatasets <= 10)
    {
        std::string dataSetName = "STFT" + std::to_string(numDatasets);
        H5::DataSpace dataspace(2, dims);
        H5::DataSet dataset = cir.createDataSet(dataSetName, H5::PredType::NATIVE_DOUBLE, dataspace);
        dataset.write(signal.stftAbs.data(), H5::PredType::NATIVE_DOUBLE);
    }
}

void IQfile(stftTest &sig, std::string fileName)
{
    std::filesystem::path backOneDir = std::filesystem::current_path().parent_path();
    std::string hdfName = fileName;
    bool test = true;
    std::string folder = test == true ? "testOutput/" : "database/";
    std::string directory = (backOneDir / folder).string();
    std::string fullPath = directory + hdfName + ".h5";
    hsize_t dims[2] = {sig.output.size(), 2};           // Change the second dimension to 2
    std::vector<double> storage(sig.output.size() * 2); // Adjust the size of the storage vector
    for (size_t i = 0; i < sig.output.size(); ++i)
    {
        storage[i * 2] = std::real(sig.ifft[i]);     // Store the real part at even indices
        storage[i * 2 + 1] = std::imag(sig.ifft[i]); // Store the imaginary part at odd indices
    }

    H5::H5File file;
    if (!std::filesystem::exists(fullPath))
    {
        file = H5::H5File(fullPath, H5F_ACC_TRUNC);
        H5::Group group = file.createGroup("/fftOutput");
    }
    else
    {
        file = H5::H5File(fullPath, H5F_ACC_RDWR);
    }
    // Open the group and count the number of datasets
    H5::Group cir = file.openGroup("/fftOutput");
    int numDatasets = cir.getNumObjs();

    if (numDatasets <= 0)
    {
        std::string dataSetName = "STFT0";
        H5::DataSpace dataspace(2, dims);
        H5::DataSet dataset = cir.createDataSet(dataSetName, H5::PredType::NATIVE_DOUBLE, dataspace);
        dataset.write(storage.data(), H5::PredType::NATIVE_DOUBLE);
    }
    else if (numDatasets >= 1 && numDatasets <= 10)
    {
        std::string dataSetName = "STFT" + std::to_string(numDatasets);
        H5::DataSpace dataspace(2, dims);
        H5::DataSet dataset = cir.createDataSet(dataSetName, H5::PredType::NATIVE_DOUBLE, dataspace);
        dataset.write(storage.data(), H5::PredType::NATIVE_DOUBLE);
    }
}

void sfft(std::vector<double> &windows, stftTest &sig, std::vector<std::complex<double>> iqBuff)
{

    using std::chrono::duration;
    using std::chrono::duration_cast;
    using std::chrono::high_resolution_clock;
    using std::chrono::milliseconds;
    auto t1 = high_resolution_clock::now();

    double signalLength = iqBuff.size();
    double windowSize = 256;
    double hopSize = windowSize / 2;
    sig.windowSize = windowSize;
    sig.hopSize = hopSize;

    double N = windowSize;
    double overlapLength = windowSize - hopSize;
    fftw_complex *in = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * windowSize);
    fftw_complex *out = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * windowSize);
    fftw_complex *origData = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * windowSize);
    fftw_plan p;
    fftw_plan ifft;
    double windowCount = floor((signalLength - windowSize) / hopSize);
    sig.windowCount = windowCount;
    sig.signalLength = signalLength;
    p = fftw_plan_dft_1d(windowSize, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    ifft = fftw_plan_dft_1d(windowSize, out, origData, FFTW_BACKWARD, FFTW_ESTIMATE);

    std::cout << " WINDOW COUNT: " << sig.windowCount << "\n";
    std::cout << " WINDOW SIZE: " << windowSize << "\n";
    std::cout << " HOP SIZE: " << hopSize << "\n";
    std::cout << "Signal Length: " << iqBuff.size() << "\n";

    hamming(windowSize, windows);
    // Process each chunk of the signal
    int chunkPosition = 0;

    int readIndex;

    // Should we stop reading in chunks?
    int bStop = 0;
    int numChunks = 0;
    std::cout << "Starting STFT\n";
    sig.output.resize(signalLength);
    sig.ifft.resize(signalLength);
    while (chunkPosition < signalLength && !bStop)
    {

        // Copy the chunk into our buffer
        for (int i = 0; i < windowSize; i++)
        {

            readIndex = chunkPosition + i;

            if (readIndex < signalLength)
            {

                // Note the windowing!
                in[i][0] = std::real(iqBuff[readIndex] * windows[i]);
                in[i][1] = std::imag(iqBuff[readIndex] * windows[i]);
            }
            else
            {

                // we have read beyond the signal, so zero-pad it!

                in[i][0] = 0.0;
                in[i][1] = 0.0;

                bStop = 1;
            }
        }

        // Perform the FFT on our chunk
        fftw_execute(p);
        for (int j = 0; j < windowSize; j++)
        {
            out[j][0] *= 1. / windowSize;
            out[j][1] *= 1. / windowSize;
        }
        fftw_execute(ifft);

        // Uncomment to see the raw-data output from the FFT calculation
        // std::cout << "Column: " << chunkPosition << std::endl;
        /*for (int i = 0; i < windowSize; i++)
        {
            fprintf(stdout, "fft_result[%d] = { %2.2f, %2.2f }\n",
                    i, out[i][0], out[i][1]);
        }*/

        // Copy the first (windowSize/2 + 1) data points into your spectrogram.
        // We do this because the FFT output is mirrored about the nyquist
        // frequency, so the second half of the data is redundant. This is how

        // Matlab's spectrogram routine works.

        for (int j = 0; j < windowSize; j++)
        {
            int outputIndex = chunkPosition + j;
            if (outputIndex >= 0 && outputIndex < signalLength)
            {
                sig.output[outputIndex] = std::complex<double>(out[j][0], out[j][1]); // Assuming real part is what we need
                sig.ifft[outputIndex] = std::complex<double>(origData[j][0] / windows[j], origData[j][1] / windows[j]);
                //  std::cout << "fft_result[" << j << "] = { " << sig.output[outputIndex] << " }\n";
            }
        }

        chunkPosition += hopSize;
        numChunks++;

    } // Excuse the formatting, the while ends here.

    fftw_destroy_plan(p);

    fftw_free(in);

    fftw_free(out);
    auto t2 = high_resolution_clock::now();

    duration<double, std::milli> ms_double = t2 - t1;

    std::cout << "Time to execute FFT for 1 thread" << ": " << ms_double.count() << "ms \n";

    std::cout << "size of output " << sig.output.size() << "\n";

    readIndex = 0;
    for (int i = 0; i < windowCount; i++)
    {
        std::vector<double> temp;
        readIndex = i * hopSize;
        for (int j = 0; j < windowSize; j++)
        {
            double tempV = std::abs(sig.output[readIndex + j]);
            temp.push_back(tempV);
            // check is abs is negative
            if (tempV < 0)
            {
                std::cout << "Negative Value: " << tempV << "\n";
            }
        }
        sig.stftAbs.push_back(temp);
    }

    std::cout << "size of stftAbs " << sig.stftAbs[0].size() << "\n";
    std::cout << "windows of stftAbs " << sig.stftAbs.size() << "\n";
}

void stftOverlapAdd(std::vector<double> &windows, stftTest &sig, std::vector<std::complex<double>> iq)
{
    sig.windowSize = 128;
    sig.hopSize = 64;
    double windowSize = sig.windowSize;
    double hopSize = sig.hopSize;
    sig.windowCount = floor((iq.size()) / hopSize);
    double signalLength = iq.size();
    std::cout << "WINDOW COUNT" << sig.windowCount << "\n";
    std::cout << "SIGNAL LENGTH" << signalLength << "\n";
    fftw_complex *in = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * windowSize);
    fftw_complex *out = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * windowSize);
    fftw_complex *ifft = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * windowSize);

    fftw_plan p;

    p = fftw_plan_dft_1d(windowSize, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

    hamming(windowSize, windows);
    sig.output.resize(windowSize * sig.windowCount);
    std::cout << "outputSize: " << sig.output.size() << "\n";
    for (int i = 0; i < signalLength; i += hopSize)
    {
        int readIndex = i;
        std::cout << "Read Index: " << readIndex << "\n";
        for (int j = 0; j < windowSize; j++)
        {
            in[j][0] = std::real(iq[readIndex + j]) * windows[j];
            in[j][1] = std::imag(iq[readIndex + j]) * windows[j];
        }

        fftw_execute(p);

        for (int j = 0; j < windowSize; j++)
        {
            if (i < signalLength)
            {

                sig.output[i + j] = std::complex<double>(out[j][0], out[j][1]);
                // std::cout << "fft_result[" << j << "] = { " << sig.output[j] << " }\n";
            }
        }
    }

    fftw_destroy_plan(p);
    fftw_free(ifft);
    fftw_free(out);
    fftw_free(in);
}

test readHdf(std::string fileName)
{

    try
    {
        test wave;
    std::filesystem::path backOneDir = std::filesystem::current_path().parent_path();
    std::string hdfName = fileName;
    bool test = true;
    std::string folder = test == true ? "testOutput/" : "testQAM/";
    std::string directory = (backOneDir / folder).string();
    std::string fullPath = directory + hdfName + ".h5";
    std::cout << "Full path to HDF5 file: " << fullPath << std::endl;

    H5::H5File file(fullPath, H5F_ACC_RDONLY);
    H5::Group group = file.openGroup("/fftOutput");
    int numDatasets = group.getNumObjs();
    std::cout << "Number of datasets: " << numDatasets << std::endl;

    if (numDatasets <= 0)
    {
        std::cout << "No datasets found in the HDF file." << std::endl;
        return wave;
    }

    std::string datasetName = "STFT0";
    if (!H5Lexists(group.getId(), datasetName.c_str(), H5P_DEFAULT))
    {
        std::cout << "Dataset " << datasetName << " does not exist." << std::endl;
        return wave;
    }

    H5::DataSet dataset = group.openDataSet(datasetName);
    H5::DataSpace dataspace = dataset.getSpace();
    int rank = dataspace.getSimpleExtentNdims();
    std::cout << "Rank of dataspace: " << rank << std::endl;
std::vector<std::vector<double>> buffer;
    hsize_t dims[2];
    dataspace.getSimpleExtentDims(dims, NULL);
    int windowSize = dims[0];
    int windowCount = dims[1];
    std::cout << "Window size: " << windowSize << ", Window count: " << windowCount << std::endl;

    std::vector<double> data(dims[0] * dims[1]);
    H5::DataSpace memspace(1, dims);
    dataset.read(data.data(), H5::PredType::NATIVE_DOUBLE);
    buffer.push_back(data);

    for(int  i = 0; i < windowCount; i++)
    {
        std::vector<double> temp;
        for(int j = 0; j < windowSize; j++)
        {
            temp.push_back(data[i * windowSize + j]);
        }
        wave.hdfAbs.push_back(temp);
    }

    
        return wave;
    } catch (H5::Exception &error)
    {
        error.printErrorStack();
    }
    catch (std::exception &e)
{
    std::cerr << "Standard exception: " << e.what() << std::endl;
}


}
