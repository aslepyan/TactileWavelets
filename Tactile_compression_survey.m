%% Load wavelets, dataset, and make range
clear all
close all
load('all_wavelets.mat')
% access wavelet by all_wavelets{i} 
% -- remove 11-45 and 60-94 and 109-114 (db10-45 and sym40-45 and mb12.3-32.3) --> too big filters
all_wavelets([11:45,60:94,109:114]) = [];

% Load STAG dataset
load('Stag.mat')
stag = double(permute(stag(:,:,:),[3,2,1])); %flop order of stuff
stag = stag(:,:,1:33550); %335450
% define filter
filt_stag = zeros(size(stag));

lpFilt = designfilt('lowpassiir','FilterOrder',5, ...
         'PassbandFrequency',1, 'SampleRate',7);

for i=1:32
    for j=1:32
        z = squeeze(stag(i,j,:));
        mz = mean(z);
        z = z-mz;
        q = [fliplr(z);z;fliplr(z)];
        y = filtfilt(lpFilt,q);
        y = y(length(z)+1:2*length(z));
        y = y+mz;
        filt_stag(i,j,:) = y;
    end
end

%plot(squeeze(stag(1,1,:)))
%hold on
%plot(squeeze(filt_stag(1,1,:)))
tacData = filt_stag;

zero_zone = zeros(32,32);
zero_zone(4:6,1:14) = 1;
zero_zone(10,1:14) = 1;
zero_zone(14:15,1:14) = 1;
zero_zone(19:32,1:25) = 1;
zero_zone(19:32,30:32) = 1;
zero_zone = zero_zone';
num_zero_t2axels = sum(zero_zone,'all'); % 476 sensors
ind_nnz = find(~zero_zone); %548 indices

clear stag i j lpFilt mz q y z

% Make range
close all
x = linspace(2,3,4);
y=x.^x;
in_min = min(y);
in_max = max(y);
out_max = 0.01;
out_min = 0.0001;
range = (y - in_min) * (out_max - out_min) / (in_max - in_min) + out_min;
range = fliplr(range);
range = [range, 0];
clear x y in_min in_max out_max out_min
%plot(range)
%% 1D
L = length(all_wavelets);
R = length(range);

j = 0;
parfor inc=1:length(range) %.01 to .0001 and 0
    mse = range(inc);
    j = j+1;
    % dimension 1
    % loop through the wavelets
    mse
    for i=1:length(all_wavelets)
        answers = cell(11,2);
        wname = all_wavelets{i};
        [sparse_rep,sparse_recon,book_keeping,quant,mses,quants, ...
            sparsity,q_max,num_nnz_bits,energy_ratio,means] = ...
            sparsify_W1_mse(tacData,mse,wname);
        % save num_nnz_bits,q_max,sparsity
        answers{1,1} = sparse_rep;
        answers{2,1} = sparse_recon;
        answers{3,1} = book_keeping;
        answers{4,1} = quant;
        answers{5,1} = mses;
        answers{6,1} = quants;
        answers{7,1} = sparsity;
        answers{8,1} = q_max;
        answers{9,1} = num_nnz_bits;
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
        answers{9,2} = 'num_nnz_bits';
        answers{10,2} = 'energy_ratio';
        answers{11,2} = 'means';

        parsave(['outputs-1d-filt/ratio_',num2str(mse),'_wavelet_',all_wavelets{i},'.mat'],answers);
        i/length(all_wavelets)
    end
end
%% 2D
L = length(all_wavelets);
R = length(range);

j = 0;
parfor inc=1:length(range) %.01 to .0001 and 0
    mse = range(inc);
    j = j+1;
    % dimension 1
    % loop through the wavelets
    mse
    for i=1:length(all_wavelets)
        answers = cell(11,2);
        wname = all_wavelets{i};
        
        [sparse_rep,sparse_recon,book_keeping,quant,mses,quants, ...
            sparsity,q_max,num_nnz_bits,energy_ratio,means] = ...
            sparsify_W2_mse(tacData,mse,wname,zero_zone);

        % save num_nnz_bits,q_max,sparsity
        answers{1,1} = sparse_rep;
        answers{2,1} = sparse_recon;
        answers{3,1} = book_keeping;
        answers{4,1} = quant;
        answers{5,1} = mses;
        answers{6,1} = quants;
        answers{7,1} = sparsity;
        answers{8,1} = q_max;
        answers{9,1} = num_nnz_bits;
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
        answers{9,2} = 'num_nnz_bits';
        answers{10,2} = 'energy_ratio';
        answers{11,2} = 'means';

        parsave(['outputs-2d-filt/ratio_',num2str(mse),'_wavelet_',all_wavelets{i},'.mat'],answers);

        i/length(all_wavelets)
    end
end
%% 3D

j = 0;
parfor inc=1:length(range) %.01 to .0001 and 0
    mse = range(inc);
    % dimension 1
    % loop through the wavelets
    mse
    for i=1:length(all_wavelets)
        answers = cell(10,2);
        wname = all_wavelets{i};

        [sparse_rep,sparse_recon,quant,mses,quants,sparsity,q_max, ...
            num_nnz_bits,energy_ratio,means] = ...
            sparsify_W3_mse(tacData,mse,wname,zero_zone);

        % save num_nnz_bits,q_max,sparsity
        answers{1,1} = sparse_rep;
        answers{2,1} = sparse_recon;
        answers{3,1} = quant;
        answers{4,1} = mses;
        answers{5,1} = quants;
        answers{6,1} = sparsity;
        answers{7,1} = q_max;
        answers{8,1} = num_nnz_bits;
        answers{9,1} = energy_ratio;
        answers{10,1} = means;

        answers{1,2} = 'sparse_rep';
        answers{2,2} = 'sparse_recon';
        answers{3,2} = 'quant';
        answers{4,2} = 'mses';
        answers{5,2} = 'quants';
        answers{6,2} = 'sparsity';
        answers{7,2} = 'q_max';
        answers{8,2} = 'num_nnz_bits';
        answers{9,2} = 'energy_ratio';
        answers{10,2} = 'means';

        parsave(['outputs-3d-filt/ratio_',num2str(mse),'_wavelet_',all_wavelets{i},'.mat'],answers);

        i/length(all_wavelets)
    end
end