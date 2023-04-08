function [sparse_rep,sparse_recon,book_keeping,quant,mses,quants,sparsity,q_max,num_nnz_bits,energy_ratio,means] = sparsify_W1_mse(rawData,mse,wname)
% This function sparsifies the rawData looking at each sensor over time
% individually (1D). It sparsifies by taking the wavelet transform and then 
% quantizing.

% calculate MSE and return energy ratio..

% INPUTS:
    % rawData = data to be sparsified; - [rows x col x time]
    % mse = acceptable error
    % wname = wavelet used; typical = 'sym3'
% OUTPUTS:
    % sparse_rep = sparse representation after quantization - [rows x col x time]
    % sparse_recon = reconstruction using sparse_rep
    % book_keeping = book keeping vector for wavelet x-form
    % quant = quantiziation value --> divide by this number and round
    % mses = evolution of mse during search
    % quants = evolution of quants during search
    % sparsity = sparsity (nnz / original size)
    % q_max = max value after quantization
    % num_nnz_bits = # of bits to represent q_max * # nnz values
    % energy_ratio = energy ratio between quant and o.g. wavelet x-form
    % means = means for each sensor

%% Run this cell and uncomment below for example
% initialize
raw = rawData;
[r,c,t] = size(rawData);
temp_recon = zeros(size(rawData));
sparse_recon = temp_recon;
mses = [];
quants = [];

if mse == 0
    tolerance = 1;
else
    tolerance = strsplit(num2str(mse),'.');
    tolerance = numel(tolerance{2});
    tolerance = 0.5*10^-(tolerance+1);
end %tolerance is one extra decimal point past mse

% De-mean rawData --> this makes energy conservation more informative
means = zeros(r,c);
for row = 1:r
    for col = 1:c
        means(row,col) = mean(rawData(row,col,:));
        rawData(row,col,:) = rawData(row,col,:) - mean(rawData(row,col,:));
    end
end

% Take transform
% rawData = [rows x col x time];
[wc,book_keeping] = wavedec(squeeze(rawData(1,1,:)),fix(log2(t)),wname);
% do the first one to get the size
new_sz = length(wc);
sparse_rep = zeros(r,c,new_sz);
for row = 1:r
    for col = 1:c
        x = squeeze(rawData(row,col,:));
        [wc,book_keeping] = wavedec(x,fix(log2(t)),wname);
        sparse_rep(row,col,:) = wc;
    end
end

original_energy = norm(sparse_rep(:)); % for energy_ratio

% Start values for binary search
mx = max(sparse_rep,[],'all'); 
mn = 1e-20;
quant = mx;

% Search for best quantizer
searching = 1;

if mse == 0 %if we want to have no error (PR)
    searching = 0;
    temp_sp = sparse_rep;
    for row = 1:r
        for col = 1:c
            x = squeeze(temp_sp(row,col,:));
            x_hat = waverec(x,book_keeping,wname);
            % add back means
            sparse_recon(row,col,:) = x_hat + means(row,col);
        end
    end
end

while (searching)
    % Quantize
    temp_sp = quant*(fix(sparse_rep./quant));

    % Reconstruct and calculate MSE
    for row = 1:r
        for col = 1:c
            x = squeeze(temp_sp(row,col,:));
            x_hat = waverec(x,book_keeping,wname);
            % add back means
            sparse_recon(row,col,:) = x_hat + means(row,col);
        end
    end

    MSE = mean(((raw(:)-sparse_recon(:)).^2)./mean(raw(:)));

    mses = [mses, MSE];
    quants = [quants, quant];

    % Update quant through binary search
    old_quant = quant;
    if MSE < mse % not enough compression, quantize more
        mn = quant;
        quant = (quant + mx)/2;
    end
    if MSE > mse % too much compression, quantize less
        mx = quant;
        quant = (quant + mn)/2;
    end

    if abs(quant - old_quant) < tolerance % if values stop changing
        searching = 0;
    end

    if (mse - tolerance < MSE) && (MSE < mse + tolerance) % if correct value is reached
        searching = 0;
    end

end

% Estimate sparsity
sparse_rep = temp_sp;
sparsity = nnz(sparse_rep) / (r*c*t);
% calculate number of nnz coefficients vs original size

% what is the maximum quantized value?
q_max = round(max(sparse_rep./quant,[],'all'));
num_nnz_bits = nnz(sparse_rep) * ceil((log2(q_max))+1);

% energy conservation
comp_power = norm(sparse_rep(:));
energy_ratio = comp_power / original_energy; %compressed power over original power
end