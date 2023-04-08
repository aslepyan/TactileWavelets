# Sparsifying Tactile Data by Wavelet Transform
This repository holds the MATLAB and python code used to generate the results in the paper: “Sparsifying Tactile Data by Wavelet Transform” by A. Slepyan, M. Zakariaie, and N. Thakor which has been submitted to the IEEE Sensors journal.
This repository also contains supplemental figures and tables not included in the main text of the paper.

Explanation of supplemental:

Explanation of code:
- Tactile_compression_survey.m contains is the main script that computes surveys the sparsifying ability of various wavelets on tactile data.
- sparsify_W1_mse.m is the function that evaluates the 1D wavelet transform, quantizes, and saves various quantifying metrics
- sparsify_W2_mse.m is the function that evaluates the 2D wavelet transform, quantizes, and saves various quantifying metrics
- sparsify_W3_mse.m is the function that evaluates the 3D wavelet transform, quantizes, and saves various quantifying metrics
- all_wavelets.mat is the data file that contains the wavelets used in this survey

The gernal methods used in the paper is as follows:
1. Compute 75 different wavelet transforms in the 1D, 2D, and 3D cases at 5 different reconstructed MSE values for a total of 1125 candidate wavelet transformers
2. Quantify each transformer with 3 qualifiers: sparsity, number of bits per pixel, energy-ratio
3. Do the inverse transform for each candidate compression transform and feed the reconstructed signal to obtain a fourth qualifier, classificaion accuracy

The STAG dataset can be found at the following link. For direct use in the above code, please extract the zipped folder, open the 'metadata.mat' file, and re-save the 'pressure' variable as a new variable 'Stag.mat'. 'Stag.mat' is used for the sparsifying survey.
