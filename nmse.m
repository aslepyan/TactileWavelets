function [nmse] = nmse(original,recon)
%original and recon should have dimensions [n x t]
nmse = mean(((original(:)-recon(:)).^2)./mean(original(:)));
end