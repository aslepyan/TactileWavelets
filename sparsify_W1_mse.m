function [sparse_rep,sparse_recon,book_keeping,quant,nmses,quants,sparsity,q_max,bpp,energy_ratio,means] = sparsify_W1_mse(Data,ind_nnz,nmse,wname)
% This function sparsifies the rawData looking at each sensor over time
% individually (1D). It sparsifies by taking the wavelet transform and then 
% quantizing.

% calculate MSE and return energy ratio..

% INPUTS:
    % rawData = data to be sparsified; - [rows x col x time]
    % nnz_indx = indexes with values
    % nmse = acceptable error
    % wname = wavelet used; example = 'sym3'
% OUTPUTS:
    % sparse_rep = sparse representation after quantization, de-meaned
    % sparse_recon = reconstruction using sparse_rep
    % book_keeping = book keeping vector for wavelet x-form
    % quant = quantiziation value --> divide by this number and round
    %                                                       towards 0
    % nmses = evolution of nmse during search
    % quants = evolution of quants during binary search
    % sparsity = sparsity (nnz / original size)
    % q_max = max value after quantization
    % bpp = average # of bits per pixel
    % energy_ratio = energy ratio between sparse_rep and original wavelet x-form
    % means = average value for each sensor

%% Run this cell and uncomment below for example
% initialize
raw = Data;
[r,c,t] = size(Data);
temp_recon = zeros(size(Data));
sparse_recon = temp_recon;
nmses = [];
quants = [];

if nmse == 0
    tolerance = 1;
else
    tolerance = strsplit(num2str(nmse),'.');
    tolerance = numel(tolerance{2});
    tolerance = 0.5*10^-(tolerance+1);
end %tolerance is one extra decimal point past mse

% De-mean rawData --> this makes energy conservation more informative
means = zeros(r,c);
for row = 1:r
    for col = 1:c
        means(row,col) = mean(Data(row,col,:));
        Data(row,col,:) = Data(row,col,:) - mean(Data(row,col,:));
    end
end

nnz_raw = reshape(raw, [r*c t]);
nnz_raw = nnz_raw(ind_nnz,:);

% Take transform
% rawData = [rows x col x time];
[wc,book_keeping] = wavedec(squeeze(Data(1,1,:)),fix(log2(t)),wname);
% do the first one to get the size
new_sz = length(wc);
sparse_rep = zeros(r,c,new_sz);
for row = 1:r
    for col = 1:c
        x = squeeze(Data(row,col,:));
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

if nmse == 0 %if we want to have no error (PR)
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

    % calculate nmse over the active region (nnz_indx)
    temp_recon_mse = reshape(sparse_recon, [r*c t]);
    temp_recon_mse = temp_recon_mse(ind_nnz,:);

    NMSE = abs(mean(((nnz_raw(:)-temp_recon_mse(:)).^2)./mean(nnz_raw(:))));
    %disp(['The NMSE for ', num2str(nmse), ' ', wname,' is ',num2str(NMSE)])

    nmses = [nmses, NMSE];
    quants = [quants, quant];

    % Update quant through binary search
    old_quant = quant;
    if NMSE < nmse % not enough compression, quantize more
        mn = quant;
        quant = (quant + mx)/2;
    end
    if NMSE > nmse % too much compression, quantize less
        mx = quant;
        quant = (quant + mn)/2;
    end

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

disp(['The final NMSE for ', num2str(nmse), ' ', wname,' is ',num2str(NMSE)])
disp(['The sparsity NMSE for ', num2str(nmse), ' ', wname,' is ',num2str(NMSE)])

% Estimate sparsity
sparse_rep = temp_sp;
sparsity = nnz(sparse_rep) / (length(sparse_rep(:)));
% calculate number of nnz coefficients

% what is the maximum quantized value?
q_max = round(max(sparse_rep./quant,[],'all'));
num_nnz_bits = nnz(sparse_rep) * ceil((log2(q_max))+1);
bpp = num_nnz_bits / length(sparse_rep(:));

% energy conservation
comp_power = norm(sparse_rep(:));
energy_ratio = comp_power / original_energy; %compressed power over original power
end