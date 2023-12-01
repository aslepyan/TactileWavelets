function [sparse_rep,sparse_recon,quant,nmses,quants,sparsity,q_max,bpp,energy_ratio,means] = sparsify_D2_mse(Data,nmse,ind_nnz)
% This function sparsifies the rawData looking at each sensor over time
% individually (1D). It sparsifies by taking the discrete cosine transform
% and then quantizing.

% INPUTS:
    % rawData = data to be sparsified; - [rows x col x time]
    % mse = acceptable error
    % zero_zone = zone on sensor not to consider
% OUTPUTS:
    % sparse_rep = sparse representation after quantization - [rows x col x time]
    % sparse_recon = reconstruction using sparse_rep
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
raw = Data;
[r,c,t] = size(Data);
temp_recon = zeros(size(Data));
sparse_recon = temp_recon;
nmses = [];
quants = [];

% nnz_indx = ~ind_nnz;
% z_indx = find(ind_nnz);
% raw_mse = reshape(raw,[r*c t]);
% raw_mse = raw_mse(nnz_indx,:);

if nmse == 0
    tolerance = 1;
else
    tolerance = strsplit(num2str(nmse),'.');
    tolerance = numel(tolerance{2});
    tolerance = 0.5*10^-(tolerance+1);
end %tolerance is one extra decimal point past mse

% De-mean rawData --> this makes energy conservation more informative
means = zeros(1,t);
for time = 1:t
    means(time) = mean(Data(:,:,time),'all');
    Data(:,:,time) = Data(:,:,time) - mean(Data(:,:,time),'all');
end
nnz_raw = reshape(raw, [r*c t]);
nnz_raw = nnz_raw(ind_nnz,:);

% Take transform
% rawData = [rows x col x time];
dctc = dct2(squeeze(Data(:,:,1)));
new_sz = length(dctc);
sparse_rep = zeros(new_sz,new_sz,t);
%%
for time = 1:t
    x = squeeze(Data(:,:,time));
    dctc = dct2(x);
    sparse_rep(:,:,time) = dctc;
    sparse_recon(:,:,time) = idct(dctc) + means(time);
end

original_energy = norm(sparse_rep(:)); % for energy_ratio

% Start values for binary search
mx = max(sparse_rep,[],'all');
mn = 1e-20;
quant = mx;

% Search for best quantizer
searching = 1;

% if nmse == 0 %if we want to have no error (PR)
%     searching = 0;
%     temp_sp = sparse_rep;
%     for row = 1:r
%         for col = 1:c
%             x = squeeze(temp_sp(row,col,:));
%             x_hat = idct(x);
%             %x_hat = waverec(x,book_keeping,wname);
%             % add back means
%             sparse_recon(row,col,:) = x_hat + means(row,col);
%         end
%     end
% end

%%
while (searching)
    % Quantize
    %quant
    temp_sp = quant*(fix(sparse_rep./quant));

    % Reconstruct and calculate MSE
    for time = 1:t
        x = squeeze(temp_sp(:,:,time));
        x_hat = idct2(x);
        temp_recon(:,:,time) =x_hat + means(time);
    end

    temp_recon_mse = reshape(temp_recon, [r*c t]);
    temp_recon_mse = temp_recon_mse(ind_nnz,:);
    
    NMSE = abs(mean(((nnz_raw(:)-temp_recon_mse(:)).^2)./mean(nnz_raw(:))));

    %    MSE = mean(((raw(:)-sparse_recon(:)).^2)./mean(raw(:)));

    nmses = [nmses, NMSE];
    quants = [quants, quant];

    % Update quant through binary search
    old_quant = quant;
    if NMSE < nmse % not enough compression, quantize more
        mn = quant;
        quant = 0.5*(quant + mx);
    end
    if NMSE > nmse % too much compression, quantize less
        mx = quant;
        quant = 0.5*(quant + mn);
    end

    %if abs(quant - old_quant) < tolerance % if values stop changing
    %    searching = 0;
    %end

    if (nmse - tolerance < NMSE) && (NMSE < nmse + tolerance) % if correct value is reached
        searching = 0;
    end

    % if search fails, if number of iterations > 50, stop searching
    if length(nmses) > 50 % if values stop changing
        searching = 0;
        quant = 1e10;
        NMSE = 1e10;
    end

end

disp(['The final NMSE for ', num2str(nmse),' 2D DCT is ',num2str(NMSE)])

sparse_recon = temp_recon;
%sparse_recon = reshape(temp_recon,[r*c t]);
%sparse_recon(z_indx,:) = 510;
%sparse_recon = reshape(sparse_recon, [r c t]);

% Estimate sparsity
sparse_rep = temp_sp;
sparsity = nnz(sparse_rep) / length(sparse_rep(:));
% calculate number of nnz coefficients vs original size

% what is the maximum quantized value?
q_max = round(max(sparse_rep./quant,[],'all'));
bpp = nnz(sparse_rep) * ceil((log2(q_max))+1);
bpp = bpp / length(sparse_rep(:));


% energy conservation
comp_power = norm(sparse_rep(:));
energy_ratio = comp_power / original_energy; %compressed power over original power

end