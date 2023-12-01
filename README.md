# Sparsifying Tactile Data by Wavelet Transform
This repository holds the MATLAB code, supplementary data, and functions used to generate the results and figures in the paper: “Sparsifying Tactile Data by Wavelet Transform” by A. Slepyan, M. Zakariaie, and N. Thakor which has been submitted to the IEEE TBioCAS journal.

Explanation of supplemental:
- TestResults.xlsx contains the full results of the test

Explanation of code:
- Survey.m is the main script that tests the sparsifying ability of the wavelet library on tactile interactions.
- Figure_generation.m is the main script to run to create the figures for the paper.
- sparsify_W1_mse.m is the function that evaluates the 1D wavelet transform, quantizes, and saves various quantifying metrics
- sparsify_W2_mse.m is the function that evaluates the 2D wavelet transform, quantizes, and saves various quantifying metrics
- sparsify_W3_mse.m is the function that evaluates the 3D wavelet transform, quantizes, and saves various quantifying metrics
- sparsify_D1_mse.m is the function that evaluates the 1D DCT, quantizes, and saves various quantifying metrics
- sparsify_D2_mse.m is the function that evaluates the 2D DCT, quantizes, and saves various quantifying metrics
- sparsify_D3_mse.m is the function that evaluates the 3D DCT, quantizes, and saves various quantifying metrics
- All other files are helper functions or data used in figure generation or the main test script

The tactile interaction data is the STAG dataset which can be found at the following link. For direct use in the above code, please extract the zipped folder, open the 'metadata.mat' file, and re-save the 'pressure' variable as a new variable 'Stag.mat'. 'Stag.mat' is used for the sparsifying survey.
http://stag.csail.mit.edu/datasets/classification_lite.zip
