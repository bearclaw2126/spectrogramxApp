
#include <iostream>
#include <fstream>
#include <complex>
#include <cstdlib>
#include <vector>
#include <iterator>
#include <fftw3.h>
#include <bits/stdc++.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <arpa/inet.h>
#include <netinet/in.h>
#include <math.h>
#include <algorithm>
#include <bits/stdc++.h>
#include <omp.h>
#include <chrono>
#include "hdf5.h"
#include "H5Cpp.h"
#define _USE_MATH_DEFINES
/**
 * @brief Lookup table for the window sizes used for the STFT.
 * @details All the window sizes are powers of 2 to speed up the FFT calculations.
 */
const int windowSizeLUT[10] = {64, 128, 256, 512, 1024, 2048, 4096};

/**
 * @struct sig
 * @brief Short Time Fourier Transform (STFT) structure to store the STFT values.
 * @details Performs the portion of storing the STFT values in a buffer.Along with the window size and the total number of windows processed.
 * @param stftBuff Buffer to store the complex values of the STFT
 * @param stftAbs Buffer to store the absolute values of the STFT
 * @param stftCIR Buffer to store the independent CIR values of the STFT
 * @param windowSize Size of the window used for STFT
 * @param windowCount Total number of windows
 * @param hopSize Overlap size of the windows
 */
typedef struct stft
{
    std::vector<std::vector<std::complex<double>>> stftBuff;
    /**
     * @var stftAbs
     * @brief Buffer to store the absolute values of the STFT
     * \label stftAbs
     */
    std::vector<std::vector<double>> stftAbs;
    /**
     * @var stftCIR
     * @brief Buffer to store the independent signal values of the STFT
     */
    std::vector<std::vector<double>> stftCIR;
    int windowSize;
    int windowCount;
    int hopSize;
} sig;

/**
 *
 * @section Functions Definitions
 * @defgroup Functions
 */

/**
 *
 * @ingroup Functions
 * @brief Creates the hamming window for the STFT.
 * @details The hamming window utlizes the hamming function and applies it over a set length the function applies is
 thus : \f[w[n] = 0.54 - 0.46 \cos(\frac{2\pi n}{(w-1))}\f]
 * @param w window size
 * @param  window reference to the window vector
 */
void hamming(int lo, std::vector<double> &c);

/**
 * @ingroup Functions
 * @brief This function calculates the Short Time Fourier Transform (STFT) of the input signal.
 * @details The function calculates the STFT of the input signal using the Fast Fourier Transform (FFT) algorithm.
 * The function uses the FFTW library to perform the FFT calculations.
 * Based off the formula for the stft which is
 * \f[x[m,n] = \sum_{k=0}^{m+(N-1)} x[n]w[k]e^{-j2\pi km/N}\f]
 * @param window The window vector used for the STFT
 * @param iq The input signal in the time domain
 * @return The STFT structure containing the STFT values
 */
stft stftCalc(std::vector<double> &c, std::vector<std::complex<double>> &iq);

void cirValue(stft &s);

void writeHdf5(stft &in, int round,std::string firstName);

void signalHandler(double signal);

void hdf5file(stft &sig, std::string fileName);

void hdfSTFTtest(stft &sig, std::string fileName);