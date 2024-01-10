# Compressing and Sparsifying Tactile Interactions by Wavelet Transform
This repository holds the MATLAB code, supplementary data, and functions used to generate the results and figures in the paper: “Compressing and Sparsifying Tactile Interactions by Wavelet Transform” by A. Slepyan, M. Zakariaie, T. Tran, and N. Thakor.

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

The tactile interaction data is the STAG dataset which can be found at the following link.
http://stag.csail.mit.edu/datasets/classification_lite.zip

Three additional (large) files are needed to run the test code and figure generation code and are found at this onedrive link, they should be placed in the same directory as the functions when executing:
[Here](https://livejohnshopkins-my.sharepoint.com/:f:/g/personal/aslepya1_jh_edu/EolkBfklbJ9PuXpYAMATCF0B4ggfkofcBk__PjZKTxfghA?e=GpTKU5)
