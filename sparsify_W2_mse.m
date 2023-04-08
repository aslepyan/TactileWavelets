function [sparse_rep,sparse_recon,book_keeping,quant,mses,quants,sparsity,q_max,num_nnz_bits,energy_ratio,means] = sparsify_W2_mse(rawData,mse,wname,zero_zone)
% This function sparsifies the rawData looking at the entire sensor array
% at each time point individually (2D). It sparsifies by taking the wavelet 
% transform and then quantizing.

% INPUTS:
    % rawData = data to be sparsified; - [rows x col x time]
    % power_ratio = acceptable power ratio; typical value is 0.999
    % wname = wavelet used; typical = 'sym3'
% OUTPUTS:
    % sparse_rep = sparse representation after quantization - [rows*col x time]
    % sparse_recon = reconstruction using sparse_rep
    % book_keeping = book keeping vector for wavelet x-form
    % quant = quantiziation value --> divide by this number and round
    % energies = evolution of energy ratio during search
    % quants = evolution of quants during search
    % sparsity = sparsity (nnz / original size)
    % q_max = max value after quantization
    % num_nnz_bits = # of bits to represent q_max * # nnz values
    % MSE = MSE between sparse_recon and and rawData
    % means = means for each frame of the sensor array over time
    
%% Run this cell and uncomment below for example

% initialize
raw = rawData;
[r,c,t] = size(rawData);
temp_recon = zeros(size(rawData));
sparse_recon = temp_recon;
mses = [];
quants = [];
%tolerance = (1-power_ratio)/10;

nnz_indx = ~zero_zone;
z_indx = find(zero_zone);
raw_mse = reshape(raw,[r*c t]);
raw_mse = raw_mse(nnz_indx,:);

if mse ==0
    tolerance = 1;
else
    tolerance = strsplit(num2str(mse),'.');
    tolerance = numel(tolerance{2});
    tolerance = 0.5*10^-(tolerance+1);
end

% De-mean rawData
means = zeros(1,t);
for time = 1:t
    means(time) = mean(rawData(:,:,time),'all');
    rawData(:,:,time) = rawData(:,:,time) - mean(rawData(:,:,time),'all');
end

% Take transform
%rawData = [rows x col x time];
[wc,book_keeping] = wavedec2(squeeze(rawData(:,:,1)),wmaxlev([r c],wname),wname);
new_sz = length(wc);
sparse_rep = zeros(new_sz,t);
for time = 1:t
    x = squeeze(rawData(:,:,time));
    [wc,book_keeping] = wavedec2(x,wmaxlev([r c],wname),wname);
    sparse_rep(:,time) = wc;
    sparse_recon(:,:,time) = waverec2(wc,book_keeping,wname) + means(time);
end

original_energy = norm(sparse_rep(:));

mx = max(sparse_rep,[],'all');
mn = 1e-20;
quant = mn;

searching = 1;

if mse == 0
    searching = 0;
    temp_sp = sparse_rep;
end

% Search for best quantizer
while (searching)
    % Quantize
    temp_sp = quant*(fix(sparse_rep./quant));

    % Reconstruct
    for time = 1:t
        x = squeeze(temp_sp(:,time));
        x_hat = waverec2(x,book_keeping,wname);
        temp_recon(:,:,time) = x_hat + means(time);
    end

    % Zero out the zero-zone
    temp_recon_mse = reshape(temp_recon, [r*c t]);
    temp_recon_mse = temp_recon_mse(nnz_indx,:);
    MSE = mean(((raw_mse(:)-temp_recon_mse(:)).^2)./mean(raw_mse(:)));
   
    mses = [mses, MSE];
    quants = [quants, quant];

    % Update quant through binary search
    old_quant = quant;
    if MSE < mse % not enough compression, quantize more
        mn = quant;
        %quant = ((iter+1)/iter)*(quant + mx)/2;
        quant = (quant + mx)/2;
        %prev_quant - (prev_quant + quant)/2;
    end
    if MSE > mse % too much compression, quantize less
        mx = quant;
        %quant = (iter/(iter+1))*(quant + mn)/2;
        quant = (quant + mn)/2;
        %prev_quant + (prev_quant + quant)/2;
    end

    if abs(quant - old_quant) < tolerance % if values stop changing
        searching = 0;
    end

    if (mse - tolerance < MSE) && (MSE < mse + tolerance) % if correct value is reached
        searching = 0;
    end
end

%set zero_zone values to 510
temp_recon = reshape(temp_recon,[r*c t]);
temp_recon(z_indx,:) = 510;
sparse_recon = reshape(temp_recon, [r c t]);

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