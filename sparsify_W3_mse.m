function [sparse_rep,sparse_recon,quant,mses,quants,sparsity,q_max,num_nnz_bits,energy_ratio,means] = sparsify_W3_mse(rawData,mse,wname,zero_zone)
% This function sparsifies the rawData looking at the entire sensor array
% over time, collectively (3D). It sparsifies by taking the wavelet 
% transform and then quantizing.

% INPUTS:
    % rawData = data to be sparsified; - [rows x col x time]
    % mse = acceptable error
    % wname = wavelet used; typical = 'sym3'
    % zero_zone = zone to not include when calculating MSE
% OUTPUTS:
    % sparse_rep = sparse representation after quantization - 1x1 struct
    %                               sparse_rep.dec{i} has decomp data
    % sparse_recon = reconstruction using sparse_rep
    % quant = quantiziation value --> divide by this number and round
    % mses = evolution of mse during search
    % quants = evolution of quants during search
    % sparsity = sparsity (nnz / original size)
    % q_max = max value after quantization
    % num_nnz_bits = # of bits to represent q_max * # nnz values
    % energy_ratio =  energy ratio between quant and o.g. wavelet x-form
    % means = mean of entire dataset
    
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

if mse==0
    tolerance = 1;
else
    tolerance = strsplit(num2str(mse),'.');
    tolerance = numel(tolerance{2});
    tolerance = 0.5*10^-(tolerance+1);
end

% De-mean rawData
means = mean(rawData,'all');
rawData = rawData - means;

mses = [];
quants = [];

% Take transform
lvl = wmaxlev(min([r,c,t]),wname);
if lvl == 0
    lvl = 1;
end
wc = wavedec3(rawData,lvl,wname);

% Find max val
mx = max(wc.dec{1},[],'all');
all_coefs = [];
for i=1:length(wc.dec)
    dec =  wc.dec{i};
    all_coefs = [all_coefs; dec(:)];
    if max(dec,[],'all') > mx
        mx = max(dec,[],'all');
    end
end

original_energy = norm(all_coefs(:));

% Reconstruct
temp_sp = wc;
temp_recon = waverec3(temp_sp);

mn = 1e-20;
quant = mn;

searching = 1;

q_max = 1;
if mse == 0
    searching = 0;
end

% Search for best quantizer
while (searching)
    % Quantize
    q_max = 0;
    all_comp_coefs = [];
    for i=1:length(temp_sp.dec)
        dec = quant*(fix(wc.dec{i}./quant));
        temp_sp.dec{i} = dec;
        all_comp_coefs = [all_comp_coefs; dec(:)];
        rr = max(temp_sp.dec{i}./quant,[],'all');
        q_max = max([q_max,rr]);
    end

    % Reconstruct
    temp_recon = waverec3(temp_sp) + means;

    % Dont inlcude zero-zone in MSE calculation
    temp_recon_mse = reshape(temp_recon, [r*c t]);
    temp_recon_mse = temp_recon_mse(nnz_indx,:);
    % Calculate MSE
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
sparse_rep = temp_sp;

%set zero_zone values to 510
temp_recon = reshape(temp_recon,[r*c t]);
temp_recon(z_indx,:) = 510;
sparse_recon = reshape(temp_recon, [r c t]);

% Estimate sparsity
sparsity = 0; %= nnz(sparse_rep) / (r*c*t);
for i=1:length(sparse_rep.dec)
    sparsity = sparsity + nnz(sparse_rep.dec{i});
end
num_nnz_bits = sparsity * ceil((log2(q_max)));
sparsity = sparsity / (r*c*t);

% Calculate energy conservation
if mse == 0
    comp_power = 1;
else
    comp_power = norm(all_comp_coefs(:));
end

energy_ratio = comp_power / original_energy; %compressed power over original power
end