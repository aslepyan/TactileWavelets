# Sparsifying Tactile Data by Wavelet Transform
This repository holds the MATLAB and python code used to generate the results in the paper: “Sparsifying Tactile Data by Wavelet Transform” by A. Slepyan, M. Zakariaie, and N. Thakor which has been submitted to the IEEE Sensors journal.
This repository also contains supplemental figures and tables not included in the main text of the paper.

Explanation of supplemental:

Explanation of code:
The code is broken up into 

-- took stag data --> computed 75 wavelet transforms and quantiation --> evaluated sparsity metrics --> inverse transformed with certain rules --> fed recon data into NN classifier directly ad measured class acc
MATLAB
- main fuction which is compressing the data by wavelet transform and quantization
- associated functions that evaluate stuff
Python
- 

Stag dataset is found at --> unzip, extract metadata and convert pressure to variable "Stag.mat" 
http://stag.csail.mit.edu/datasets/classification_lite.zip
