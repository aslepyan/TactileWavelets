% This is the main code file to run for the paper:
% Sparsifying Tactile Sensor Array Data by Wavelet Transform
% The code tests a library of discrete wavelet transforsms and also the
% discrete cosine transform as 1D, 2D, and 3D transforms for sparsifying
% tactile data
% The process of the survey is as follows:
% (1) Load sensor array data
% (2) Specify desired NMSE values
% (3) Transform data using candidate transform, and quantize until
%     reconstruction error matchs desired NMSE
% (4) Save reconstructed signal and calculate quantifying metrics:
%       Sparsity, Average Bits Per Pixel, Energy Ratio
%% Load data, load wavelets, and create range of NMSEs
clear all
clc
close all
load('all_wavelets.mat')
% access wavelet by all_wavelets{i}
% -- remove 11-45 and 60-94 and 109-114 (db10-45 and sym40-45 and mb12.3-32.3) --> big filters
all_wavelets([11:45,60:94,109:114]) = [];
load('filt_stag.mat')
load('zero_zone.mat') % the STAG dataset contains pixels that are not active,
% we refer to these as the 'zero_zone'
ind_nnz = find(~zero_zone); %548 indices
range = [0.01, 0.0043, 0.0015, 0.0001]; %range of NMSE values to test
%% What values are we trying to save for each run?
% vary dimension, and wavelet, and nmse
% save: sparsity, bpp, ER, sparse reps and recons
%% 1D survey

L = length(all_wavelets);
R = length(range);

j = 0;
%summary of wavelet results
wSs = zeros(length(range),length(all_wavelets));
wBPPs = zeros(length(range),length(all_wavelets));
wERs = zeros(length(range),length(all_wavelets));
wNames = cell(length(range),length(all_wavelets));

%summary of dct results
dSs = zeros(length(range),1);
dBPPs = zeros(length(range),1);
dERs = zeros(length(range),1);
dNames = cell(length(range),1);

parfor inc=1:length(range) %.01 to .0001 and 0
    nmse = range(inc);
    j = j+1;
    % dimension 1
    % loop through the wavelets
    disp(['The current nmse is ',num2str(nmse)]) %print to see progress
    for i=1:L
        answers = cell(11,2);
        wname = all_wavelets{i};
        [sparse_rep,sparse_recon,book_keeping,quant,nmses,quants, ...
            sparsity,q_max,bpp,energy_ratio,means] = ...
            sparsify_W1_mse(filt_stag,ind_nnz,nmse,wname);

        % save num_nnz_bits,q_max,sparsity
        wSs(inc,i) = sparsity;
        wBPPs(inc,i) = bpp;
        wERs(inc,i) = energy_ratio;
        wNames{inc,i} = wname + "1D" + num2str(nmse);

        answers{1,1} = sparse_rep;
        answers{2,1} = sparse_recon;
        answers{3,1} = book_keeping;
        answers{4,1} = quant;
        answers{5,1} = nmses;
        answers{6,1} = quants;
        answers{7,1} = sparsity;
        answers{8,1} = q_max;
        answers{9,1} = bpp;
        answers{10,1} = energy_ratio;
        answers{11,1} = means;

        answers{1,2} = 'sparse_rep';
        answers{2,2} = 'sparse_recon';
        answers{3,2} = 'book_keeping';
        answers{4,2} = 'quant';
        answers{5,2} = 'nmses';
        answers{6,2} = 'quants';
        answers{7,2} = 'sparsity';
        answers{8,2} = 'q_max';
        answers{9,2} = 'bpp';
        answers{10,2} = 'energy_ratio';
        answers{11,2} = 'means';

        %parsave(['outputs_1D/nmse_',num2str(nmse),'_wavelet_',all_wavelets{i},'.mat'],answers);

    end %finish wavelet loop

    % start dcts
    answers = cell(10,2);
    [sparse_rep,sparse_recon,quant,nmses,quants,sparsity,q_max,bpp, ...
        energy_ratio,means] = ...
        sparsify_D1_mse(filt_stag,nmse,ind_nnz);

    % save num_nnz_bits,q_max,sparsity
    dSs(inc) = sparsity;
    dBPPs(inc) = bpp;
    dERs(inc) = energy_ratio;
    dNames{inc} = "DCT1D" + num2str(nmse);

    answers{1,1} = sparse_rep;
    answers{2,1} = sparse_recon;
    answers{3,1} = quant;
    answers{4,1} = nmses;
    answers{5,1} = quants;
    answers{6,1} = sparsity;
    answers{7,1} = q_max;
    answers{8,1} = bpp;
    answers{9,1} = energy_ratio;
    answers{10,1} = means;

    answers{1,2} = 'sparse_rep';
    answers{2,2} = 'sparse_recon';
    answers{3,2} = 'quant';
    answers{4,2} = 'nmses';
    answers{5,2} = 'quants';
    answers{6,2} = 'sparsity';
    answers{7,2} = 'q_max';
    answers{8,2} = 'bpp';
    answers{9,2} = 'energy_ratio';
    answers{10,2} = 'means';

    %parsave(['outputs_1D/nmse_',num2str(nmse),'_DCT_.mat'],answers);

end
Ss = [wSs,dSs];
BPPs = [wBPPs,dBPPs];
ERs = [wERs,dERs];
Names = [wNames,dNames];

%save('outputs_1D/summary.mat','Ss','BPPs','ERs','Names');
%% 2D Survey
L = length(all_wavelets);
R = length(range);

j = 0;
%summary of wavelet results
wSs = zeros(length(range),length(all_wavelets));
wBPPs = zeros(length(range),length(all_wavelets));
wERs = zeros(length(range),length(all_wavelets));
wNames = cell(length(range),length(all_wavelets));

%summary of dct results
dSs = zeros(length(range),1);
dBPPs = zeros(length(range),1);
dERs = zeros(length(range),1);
dNames = cell(length(range),1);

parfor inc=1:length(range) %.01 to .0001 and 0
   nmse = range(inc);
    j = j+1;
    % dimension 1
    % loop through the wavelets
    disp(['The current nmse is ',num2str(nmse)]) %print to see progress

    for i=1:L
        answers = cell(11,2);
        wname = all_wavelets{i};
        [sparse_rep,sparse_recon,book_keeping,quant,mses,quants, ...
            sparsity,q_max,bpp,energy_ratio,means] = ...
            sparsify_W2_mse(filt_stag,nmse,wname,zero_zone);
        
        wSs(inc,i) = sparsity;
        wBPPs(inc,i) = bpp;
        wERs(inc,i) = energy_ratio;
        wNames{inc,i} = wname + "2D" + num2str(nmse);

        answers{1,1} = sparse_rep;
        answers{2,1} = sparse_recon;
        answers{3,1} = book_keeping;
        answers{4,1} = quant;
        answers{5,1} = mses;
        answers{6,1} = quants;
        answers{7,1} = sparsity;
        answers{8,1} = q_max;
        answers{9,1} = bpp;
        answers{10,1} = energy_ratio;
        answers{11,1} = means;

        answers{1,2} = 'sparse_rep';
        answers{2,2} = 'sparse_recon';
        answers{3,2} = 'book_keeping';
        answers{4,2} = 'quant';
        answers{5,2} = 'mses';
        answers{6,2} = 'quants';
        answers{7,2} = 'sparsity';
        answers{8,2} = 'q_max';
        answers{9,2} = 'bpp';
        answers{10,2} = 'energy_ratio';
        answers{11,2} = 'means';

        %parsave(['outputs_3D/nmse_',num2str(mse),'_wavelet_',all_wavelets{i},'.mat'],answers);
    
    end %finish wavelet loop

    [sparse_rep,sparse_recon,quant,nmses,quants,sparsity,q_max,bpp, ...
        energy_ratio,means] = sparsify_D2_mse(filt_stag,nmse,ind_nnz);

    dSs(inc) = sparsity;
    dBPPs(inc) = bpp;
    dERs(inc) = energy_ratio;
    dNames{inc} = "DCT2D" + num2str(nmse);

    answers = cell(10,2);

    answers{1,1} = sparse_rep;
    answers{2,1} = sparse_recon;
    answers{3,1} = quant;
    answers{4,1} = nmses;
    answers{5,1} = quants;
    answers{6,1} = sparsity;
    answers{7,1} = q_max;
    answers{8,1} = bpp;
    answers{9,1} = energy_ratio;
    answers{10,1} = means;

    answers{1,2} = 'sparse_rep';
    answers{2,2} = 'sparse_recon';
    answers{3,2} = 'quant';
    answers{4,2} = 'nmses';
    answers{5,2} = 'quants';
    answers{6,2} = 'sparsity';
    answers{7,2} = 'q_max';
    answers{8,2} = 'bpp';
    answers{9,2} = 'energy_ratio';
    answers{10,2} = 'means';

    %parsave(['outputs_2D/nmse_',num2str(nmse),'_DCT_.mat'],answers);
end
Ss = [wSs,dSs];
BPPs = [wBPPs,dBPPs];
ERs = [wERs,dERs];
Names = [wNames,dNames];

%save('outputs_2D/summary.mat','dSs','dBPPs','dERs','dNames');


%% 3D Survey
%just do DCT for now... 

dSs = zeros(length(range),1);
dBPPs = zeros(length(range),1);
dERs = zeros(length(range),1);
dACCs = zeros(length(range),1);
dNames = cell(length(range),1);
accuracy = 0;

parfor inc=1:length(range) %.01 to .0001 and 0
        nmse = range(inc);
    j = j+1;
    % loop through the wavelets
    disp(['The current nmse is ',num2str(nmse)]) %print to see progress
    for i=1:L
        answers = cell(10,2);
        wname = all_wavelets{i};

        [sparse_rep,sparse_recon,quant,mses,quants,sparsity,q_max, ...
            bpp,energy_ratio,means] = ...
            sparsify_W3_mse(filt_stag,nmse,wname,zero_zone);

        wSs(inc,i) = sparsity;
        wBPPs(inc,i) = bpp;
        wERs(inc,i) = energy_ratio;
        wNames{inc,i} = wname + "3D" + num2str(nmse);

        answers{1,1} = sparse_rep;
        answers{2,1} = sparse_recon;
        answers{3,1} = quant;
        answers{4,1} = mses;
        answers{5,1} = quants;
        answers{6,1} = sparsity;
        answers{7,1} = q_max;
        answers{8,1} = bpp;
        answers{9,1} = energy_ratio;
        answers{10,1} = means;

        answers{1,2} = 'sparse_rep';
        answers{2,2} = 'sparse_recon';
        answers{3,2} = 'quant';
        answers{4,2} = 'mses';
        answers{5,2} = 'quants';
        answers{6,2} = 'sparsity';
        answers{7,2} = 'q_max';
        answers{8,2} = 'bpp';
        answers{9,2} = 'energy_ratio';
        answers{10,2} = 'means';

        %parsave(['outputs_3D/nmse_',num2str(nmse),'_wavelet_',all_wavelets{i},'.mat'],answers);
    
    end %finish wavelet loop
    
    %start DCTs
    [sparse_rep,sparse_recon,quant,nmses,quants,sparsity,q_max,bpp, ...
    energy_ratio,means] = sparsify_D3_mse(filt_stag,nmse,ind_nnz);

    dSs(inc) = sparsity;
    dBPPs(inc) = bpp;
    dERs(inc) = energy_ratio;
    dACCs(inc) = accuracy;
    dNames{inc} = "DCT3D" + num2str(nmse);

    answers = cell(10,2);

    answers{1,1} = sparse_rep;
    answers{2,1} = sparse_recon;
    answers{3,1} = quant;
    answers{4,1} = nmses;
    answers{5,1} = quants;
    answers{6,1} = sparsity;
    answers{7,1} = q_max;
    answers{8,1} = bpp;
    answers{9,1} = energy_ratio;
    answers{10,1} = means;

    answers{1,2} = 'sparse_rep';
    answers{2,2} = 'sparse_recon';
    answers{3,2} = 'quant';
    answers{4,2} = 'nmses';
    answers{5,2} = 'quants';
    answers{6,2} = 'sparsity';
    answers{7,2} = 'q_max';
    answers{8,2} = 'bpp';
    answers{9,2} = 'energy_ratio';
    answers{10,2} = 'means';

    %parsave(['outputs_3D/nmse_',num2str(nmse),'_DCT_.mat'],answers);
end
Ss = [wSs,dSs];
BPPs = [wBPPs,dBPPs];
ERs = [wERs,dERs];
Names = [wNames,dNames];
%save('outputs_3D/summary.mat','dSs','dBPPs','dERs','dNames');