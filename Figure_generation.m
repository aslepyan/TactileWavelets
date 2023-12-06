%% This file is used to generate the figures for the paper
% Run each section at a time to generate each figure
% Figures presented in the paper are resized/ordered as needed
%% Figure 1a
load('filt_stag.mat')
rawData = filt_stag;
clear stag
figure()
reshaped_raw = reshape(rawData, [1024 33550]);
nnz_entries = find(reshaped_raw(:,1) > 510);
region = reshaped_raw(nnz_entries,1:100);
region = 1023.*(region - min(min(region))) ./ (max(max(region)) - min(min(region)));
imagesc(region)
ylabel('Taxel Number')
xlabel('Time')
title('Applied Pressure over Time')
colorbar()
figure()
thresh_region = region;
thresh_region(thresh_region<150) = 0;
imagesc(thresh_region)
ylabel('Taxel Number')
xlabel('Time')
title('Reconstructed Pressure over Time')
colorbar()
%% Figure 1b and Fig 2
load('filt_stag.mat')

rawData = filt_stag;
load('zero_zone.mat')

close all
wname = 'haar';
img1 = filt_stag(:,:,17561) - 510;
img2 = filt_stag(:,:,17565) - 510;
img3 = filt_stag(:,:,17570) - 510;

img = img1;

[C,S] = wavedec2(img,2,wname);

[h1,v1,d1] = detcoef2('all',C,S,1);
[h2,v2,d2] = detcoef2('all',C,S,2);
A = appcoef2(C,S,wname,2);
% 
dwt_img = zeros(32,32);
dwt_img(1:8,1:8) = A;

dwt_img(9:16,1:8) = h2;
dwt_img(1:8,9:16) = v2;
dwt_img(9:16,9:16) = d2;

dwt_img(17:32,1:16) = h1;
dwt_img(1:16,17:32) = v1;
dwt_img(17:32,17:32) = d1;

dwt_img = abs(dwt_img);

subplot(1,3,3)
imagesc(dwt_img)
title('2D DWT')
% 
quant = 10;
C = quant*(fix(C./quant));

[h1,v1,d1] = detcoef2('all',C,S,1);
[h2,v2,d2] = detcoef2('all',C,S,2);
A = appcoef2(C,S,wname,2);

dwt_img = zeros(32,32);
dwt_img(1:8,1:8) = A;

dwt_img(9:16,1:8) = h2;
dwt_img(1:8,9:16) = v2;
dwt_img(9:16,9:16) = d2;

dwt_img(17:32,1:16) = h1;
dwt_img(1:16,17:32) = v1;
dwt_img(17:32,17:32) = d1;

dwt_img = abs(dwt_img);

subplot(1,3,2)
imagesc(dwt_img)
title('Quantized 2D DWT')
subplot(1,3,1)
imagesc(img)
title('Raw Tactile')

%Recon figure
%close all
xhat = waverec2(C,S,wname);
xhat(find(zero_zone)) = 0;
figure()
imagesc(xhat)
%% Fig 3 is made in powerpoint
%% Fig 4
% there are 12 families
sub_wavelets = {'db3','coif2','sym4','fk4','bl7','mb4.2','beyl','vaid','han3.3','dmey','bior2.4','rbio1.3'};
for i=1:length(sub_wavelets)
    [~,psi,xval] = wavefun(sub_wavelets{i},5);
    subplot(3,4,i)
    plot(psi)
    title(sub_wavelets{i})
    xlim([1 length(psi)])
    ylim([min(psi) max(psi)])
end
%% Fig 5
close all
clear all
load('filt_stag.mat')
w=1;
reshaped_raw = reshape(filt_stag, [1024 33550]);
%sample = reshaped_raw(261,4800:7000);
sample = reshaped_raw(261,4800:5200);
avg = mean(sample);
sample = sample - avg;
wname = 'db1';
[wc,book_keeping] = wavedec(sample,2,wname);
mse = 0.01;

subplot(4,2,1)
sig = sample;
plot(sig,'LineWidth',w)
xlim([1 length(sig)])
ylim([min(sig) max(sig)])
title('Original Signal')

subplot(4,2,3)
sig = wc(1:book_keeping(1));
plot(sig,'LineWidth',w)
xlim([1 length(sig)])
ylim([min(sig) max(sig)])
title('Approximation 2')

subplot(4,2,5)
sig = wc(book_keeping(1)+1:book_keeping(1)+book_keeping(2));
plot(sig,'LineWidth',w)
xlim([1 length(sig)])
ylim([min(sig) max(sig)])
title('Details 2')

subplot(4,2,7)
sig = wc(book_keeping(1)+book_keeping(2)+1:end);
plot(sig,'LineWidth',w)
xlim([1 length(sig)])
ylim([min(sig) max(sig)])
title('Details 1')

% Search for optimal Q

tolerance = strsplit(num2str(mse),'.');
tolerance = numel(tolerance{2});
tolerance = 0.5*10^-(tolerance+1);

mx = max(wc);
mn = 1e-20;
quant = mn;

% Search for best quantizer
searching = 1;

while (searching)
    % Quantize
    x = quant*(fix(wc./quant));

    % Reconstruct and calculate MSE
    x_hat = waverec(x,book_keeping,wname);
    % add back means

    MSE = mean(((sample-x_hat).^2))/avg;

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
temp_sp = x;
sparsi = nnz(temp_sp) / (length(sample));

q_max = round(max(temp_sp./quant,[],'all'));
num_nnz_bits = nnz(temp_sp) * ceil((log2(q_max))+1);
bpp = num_nnz_bits / length(temp_sp(:));

subplot(4,2,4)
sig = temp_sp(1:book_keeping(1));
plot(sig,'LineWidth',w)
xlim([1 length(sig)])
ylim([min(sig) max(sig)])
title('Quantized Approximation 2')

subplot(4,2,6)
sig = temp_sp(book_keeping(1)+1:book_keeping(1)+book_keeping(2));
plot(sig,'LineWidth',w)
xlim([1 length(sig)])
%ylim([min(sig) max(sig)])
title('Quantized Details 2')

subplot(4,2,8)
sig = temp_sp(book_keeping(1)+book_keeping(2)+1:end);
plot(sig,'LineWidth',w)
xlim([1 length(sig)])
%ylim([min(sig) max(sig)])
title('Quantized Details 1')

subplot(4,2,2)
plot(x_hat,'LineWidth',w)
xlim([1 length(x_hat)])
ylim([min(x_hat) max(x_hat)])
title('Sparse Reconstruction')

%top = 'Wavelet Compression by Quantization';
%subt = ['Wavelet = ',wname, '; NMSE = ',num2str(mse,2), '; Sparsity = ',num2str(sparsi,2),'; BPP = ',num2str(bpp,3)];
%sgt = sgtitle({['{\bf\fontsize{12}' top '}'],['{\fontsize{12}' subt '}']});\

%% Fig 6
clear all
load('master_data_fixed.mat')
family_names = {'beyl','bior','bl','coif','db','dmey','fk','han','mb', ...
    'rbio','sym','vaid'};

load('filt_stag.mat')
load('ind_nnz.mat')
startt = 12364; 
endt = 12412;
segment = filt_stag(:,:,startt:endt); %4900, 5030
v_segment = reshape(segment, [32*32 length(segment(1,1,:))]);
vnnz_segment = v_segment(ind_nnz,:);% - 512;

best_in_each_family = [];
for f = 1:length(family_names)
    fam = family_names{f};
    indx = [];
    spars = [];
    for i=2:900
        family = master_data{i,6};
        if isequal(family,fam)
            % save index into list
            indx = [indx, i];
            spars = [spars, master_data{i,3}];
        end
    end
    [u,v] = min(spars);
    sparsest = indx(v);
    best_in_each_family = [best_in_each_family, sparsest];
end
wavelets = cell(1,12);
for i=1:12
    wavelets{i} = [master_data{best_in_each_family(i),6},num2str(master_data{best_in_each_family(i),7})];
end

% for each wavelet in wavelets, calculate the WT of each sensor, save the
% answer into a cell (12,548)
cs = cell(12,548);
ls = cell(12,548);
for i=1:12
    wname = wavelets{i};
    for j = 1:548
        sensor = vnnz_segment(j,:);
        lvl = fix(log2(length(sensor)));
        [c,l] = wavedec(sensor,lvl,wname);
        cs{i,j} = c;
        ls{i,j} = l;
    end
end

% (a) plot MSE vs num coeffs
%max number of coefficients is determined by the length of the data =
%   length of segment = endt - startt + 1
mses = zeros(endt - startt + 1, 13);
psnrs = zeros(endt - startt + 1, 13);

% loop through each wavelet
for wav = 1:12
    for i = 1:endt - startt + 1 %loop through number to truncate
        recA = zeros(548,endt - startt + 1);
        for j = 1:548 %loop through each sensor 
            % truncate coeffs
            % mink --> find the minimum values and set those to 0 
            wt = cs{wav,j};
            K = length(wt)-i;
            [b,ind] = mink(abs(wt),K);
            wt(ind) = 0;
            recv = waverec(wt,ls{wav,j},wavelets{wav}); % reconstruct
            recA(j,:) = recv;
        end
        % calculate NMSE and PSNR over entire array
        NMSE = nmse(recA,vnnz_segment);
        PSNR = 10*log10(1023^2 / (norm(recA - vnnz_segment)^2 / (548)^2));
        % save in array --> length of data x num wavelets
        mses(i,wav) = NMSE;
        psnrs(i,wav) = PSNR;
        i
    end
    wav
end

N=15;
hex = ['#1f77b4';'#aec7e8';'#ff7f0e';'#ffbb78';'#2ca02c';'#98df8a';'#d62728';'#ff9896';'#9467bd';'#c5b0d5';'#8c564b';'#c49c94';'#e377c2';'#f7b6d2';'#7f7f7f';'#c7c7c7';'#bcbd22';'#dbdb8d';'#17becf';'#9edae5'];
raw = sscanf(hex.','#%2x%2x%2x',[3,Inf]).';

figure()

wavs = [2,5,7,8,10,11];
for i = wavs
%for i=1:12 %shape determined by family, edge color is Red (1D), face color is filter size
    wav = wavelets{i};
    semilogy(mses(:,i),'LineWidth',5)
    hold on
end
xlim([1 49])
xlabel('Number of Coefficients')
ylabel('NMSE')
%title('nmse vs num coeffs')

% (b) plot PSNR vs num coeffs
%close all
figure()
for i = wavs
%for i=1:12 %shape determined by family, edge color is Red (1D), face color is filter size
    wav = wavelets{i};
    semilogy(psnrs(:,i),'LineWidth',5)
    hold on
end
xlim([1 49])
xlabel('Number of Coefficients')
ylabel('PSNR')
%legend(wavelets,'Location','bestoutside')
%legend(wavelets(wavs),'Location','bestoutside','Orientation','horizontal','NumColumns',6)
%title('psnr vs num coeffs')

% (c) show the magnitude of sorted coefficients for each case
figure()
sorted_coeffs = cell(12,1);
for i = wavs
%for i = 1:12
    sorted_coeffs{i} = cs{i,1};
    for j = 2:548
        sorted_coeffs{i} = sorted_coeffs{i} + cs{i,j};
    end
end
 
for i = wavs
%for i=1:12
semilogy(sort(abs(sorted_coeffs{i}),'descend'),'LineWidth',5)
hold on
end
legend(wavelets(wavs),'Location','bestoutside','Orientation','horizontal','NumColumns',3)
xlim([1 49])

ylabel('Coefficient Magnitude')
xlabel('Coefficient Number')
%% Figure 7 and Fig 8p2
clear all
close all
load('master_data_fixed.mat')


family_names = {'beyl','bior','bl','coif','db','dmey','fk','han','mb', ...
    'rbio','sym','vaid'};

shape = {"+","pentagram","v","square","o","x","^","<",">","hexagram", ...
    "diamond","*"}; % family name
edgecolors = {'red','green','blue'}; %1D, 2D, 3D
%facecolors = [colormap(spring(5));colormap(turbo(5));colormap(cool(5)); ];
N=15;
hex = ['#1f77b4';'#aec7e8';'#ff7f0e';'#ffbb78';'#2ca02c';'#98df8a';'#d62728';'#ff9896';'#9467bd';'#c5b0d5';'#8c564b';'#c49c94';'#e377c2';'#f7b6d2';'#7f7f7f';'#c7c7c7';'#bcbd22';'#dbdb8d';'#17becf';'#9edae5'];
raw = sscanf(hex.','#%2x%2x%2x',[3,Inf]).';
num = size(raw,1);
facecolors = raw(1+mod(0:N-1,num),:) / 255;

MSEs = [0.01,0.0043,0.0015,0.0001];

Sparities = [0;cell2mat(master_data(2:end,3))];
nmses = [0;cell2mat(master_data(2:end,10))];

close all
margin = 5e-4*1.2;

d1ms = [];
d1ss = [];

d2ms = [];
d2ss = [];

d3ms = [];
d3ss = [];

for j = 1:4
f = figure;
nump = 0;
xs = [];
ys = [];
ds = [];
for i=2:900 %900
    if (master_data{i, 10} == MSEs(j)) %if MSE = 0.01
        %if family name = ... set the shape (6)
        x = strmatch(master_data{i, 6}, family_names);
        %if filter size = ... set the face color (8)
        %if dimension = ... set the edge color (9)
        ds = [ds, master_data{i, 9}]; % dimension
        xs = [xs, Sparities(i)];
        ys = [ys, nmses(i)];
    end
end
d1max = max(xs(ds==1));
d1min = min(xs(ds==1));
d1mean = mean(xs(ds==1));
d1ms = [d1ms, d1mean];
d1std = std(xs(ds==1));
d1ss = [d1ss, d1std];


d2max = max(xs(ds==2));
d2min = min(xs(ds==2));
d2mean = mean(xs(ds==2));
d2ms = [d2ms, d2mean];
d2std = std(xs(ds==2));
d2ss = [d2ss, d2std];

d3max = max(xs(ds==3));
d3min = min(xs(ds==3));
d3mean = mean(xs(ds==3));
d3ms = [d3ms, d3mean];
d3std = std(xs(ds==3));
d3ss = [d3ss, d3std];


s = 2;
x = [MSEs(j)-margin MSEs(j)-margin MSEs(j)+margin MSEs(j)+margin];
y = [d1mean-d1std/s d1mean+d1std/s d1mean+d1std/s d1mean-d1std/s];
y(y<0) = 1e-3;
%y = [d1min d1max d1max d1min];
hold on
fill(x,y,'red','FaceAlpha',0.3,'EdgeColor','red')
yline(d1mean,'red','LineWidth',3)


x = [MSEs(j)-margin MSEs(j)-margin MSEs(j)+margin MSEs(j)+margin];
y = [d2mean-d2std/s d2mean+d2std/s d2mean+d2std/s d2mean-d2std/s];
y(y<0) = 1e-3;
%y = [d2min d2max d2max d2min];
hold on
fill(x,y,'green','FaceAlpha',0.3,'EdgeColor','green')
yline(d2mean,'green','LineWidth',3)

x = [MSEs(j)-margin MSEs(j)-margin MSEs(j)+margin MSEs(j)+margin];
y = [d3mean-d3std/s d3mean+d3std/s d3mean+d3std/s d3mean-d3std/s];
y(y<0) = 1e-3;
%y = [d3min d3max d3max d3min];
hold on
fill(x,y,'blue','FaceAlpha',0.3,'EdgeColor','blue')
yline(d3mean,'blue','LineWidth',3)

nump = 1;
spacing = -5e-4:2e-5:5e-4;

for i=2:900 %900
    if (master_data{i, 10} == MSEs(j)) %if MSE = 0.01
        %if family name = ... set the shape (6)
        x = strmatch(master_data{i, 6}, family_names);
        %if filter size = ... set the face color (8)
        %if dimension = ... set the edge color (9)
        ds = [ds, master_data{i, 9}]; % dimension
        xs = [xs, Sparities(i)];
        ys = [ys, nmses(i)];
        plot(nmses(i) + spacing(1+mod(nump,length(spacing))), Sparities(i), shape{x}, ...
            'MarkerSize',17, 'MarkerEdgeColor', ...
            edgecolors{master_data{i, 9}}, 'MarkerFaceColor', ...
            facecolors(master_data{i, 8},:),'LineWidth', 1);
        % marker size = 20
        hold on
        nump = nump+1;
    end
end


set(gca,'yscale','log')
%set(gca,'xscale','log')
xlabel('NMSE')
ylabel('Sparsity')

margin = 5e-4*1.2;
xlim([MSEs(j)-margin MSEs(j)+margin])
%ylim([3e-3 2])

%ylim([3e-3 1.4])
%ylim([.144 .152])
f.Position = [100 100 280 700];
end

figure()
m=5;
plot(MSEs,d1ms,'red','LineWidth',m)
hold on
d1dev = d1ms-d1ss;
d1dev(d1dev<0) = 1e-3;
patch([MSEs flip(MSEs)], [d1dev flip(d1ms+d1ss)], 'red', 'FaceAlpha',0.25, 'EdgeColor','none')
plot(MSEs,d2ms,'green','LineWidth',m)
hold on
patch([MSEs flip(MSEs)], [d2ms-d2ss flip(d2ms+d2ss)], 'green', 'FaceAlpha',0.25, 'EdgeColor','none')
plot(MSEs,d3ms,'blue','LineWidth',m)
hold on
patch([MSEs flip(MSEs)], [d3ms-d3ss flip(d3ms+d3ss)], 'blue', 'FaceAlpha',0.25, 'EdgeColor','none')

set(gca,'xscale','log')
set(gca,'yscale','log')
xlabel('NMSE')
ylabel('Sparsity')
legend('1D \mu','1D \sigma','2D \mu','2D \sigma','3D \mu','3D \sigma','Location','best')
ax=gca;
k=.025;
ax.TickLength = [k, k];

%%  Fig 8 p1
clear all
close all
% best in each dimension: sym4, db1, mb4.2... pull the numbers from the
% table?
w=4;
k=.025;

nmses = [0.01, 0.0043, 0.0015, 0.0001];
sym4_1d_s = [0.00497367779531947,0.0140489180906389,0.0358596361905646,0.158693815007429];
db1_2d_s = [0.086048051, 0.157822221963487, 0.255047940340909, 0.509747985748882];
mb42_3d_s = [0.0298093156203428, 0.0724540098733234, 0.135201483327124, 0.335553686195976];

sym4_1d__bpp = [0.049736778, 0.154538099, 0.430315634, 2.380407225];
db1_2d_bpp = [1.204672713,2.525155551, 4.590862926, 11.21445569];
mb42_3d_bpp = [0.357711787, 1.159264158, 2.4336267, 7.382181096];

subplot(1,2,1)
plot(nmses,sym4_1d_s,'red','LineWidth',w)
hold on
plot(nmses,db1_2d_s,'green','LineWidth',w)
plot(nmses,mb42_3d_s,'blue','LineWidth',w)
set(gca,'xscale','log')
set(gca,'yscale','log')
xlabel('NMSE')
ylabel('Sparsity')
legend('SYM4 1D','DB1 2D','MB4.2 3D','Location','best')
ax=gca;
ax.TickLength = [k, k];

subplot(1,2,2)
plot(nmses,sym4_1d__bpp,'red','LineWidth',w)
hold on
plot(nmses,db1_2d_bpp,'green','LineWidth',w)
plot(nmses,mb42_3d_bpp,'blue','LineWidth',w)
set(gca,'xscale','log')
set(gca,'yscale','log')
xlabel('NMSE')
ylabel('Bits Per Pixel')
legend('SYM4 1D','DB1 2D','MB4.2 3D','Location','best')
ax=gca;
ax.TickLength = [k, k];



%% Fig 9
clear all
close all
load('filt_stag.mat')
load('ind_nnz.mat')
load('zero_zone.mat')

startt = 12364;
endt = 12412;

%(a)
subplot(2,4,1)

az = -37.5;
el = 15;

examp = 19/7;%25;
exampi = 20;

segment = filt_stag(:,:,startt:endt); %4900, 5030
v_segment = reshape(segment, [32*32 length(segment(1,1,:))]);
vnnz_segment = v_segment(ind_nnz,:);% - 512;
sq_seg = zeros(224,224);
num = 1;
for i=1:7
    for j=1:7
        sq_seg(32*(i-1)+1:32*(i),32*(j-1)+1:32*(j)) = segment(:,:,num);
        num = num + 1;
    end
end

t = 0:1/7:48/7;
plot(t,vnnz_segment')
xline(examp,'--r');
xlim([0 t(end)])
ylabel('ADC Reading')
xlabel('Time (sec)')
title('Original')

norm=vnnz_segment;
norm = norm - min(norm,[],'all');
norm = norm ./ (max(norm,[],'all'));
e_overall = entropy(norm);
e_temporal = 0;
for i = 1:548 %loop through each sensor
    e_temporal = e_temporal + entropy(norm(i,:));
end
e_spatial = 0;
for i = 1:49 %loop through each sensor
    e_spatial = e_spatial + entropy(norm(:,i));
end

%(b) best 1d =  SYM 4
load('nmse_0.01_wavelet_sym4_1D');
sym4_1 = answers{2,1};
sym4_1 = sym4_1(:,:,startt:endt); %4900, 5030
examp_sym4_1 = sym4_1(:,:,exampi);
sq_mb42_1 = zeros(224,224);
num = 1;
for i=1:7
    for j=1:7
        sq_mb42_1(32*(i-1)+1:32*(i),32*(j-1)+1:32*(j)) = sym4_1(:,:,num);
        num = num + 1;
    end
end
sym4_1 = reshape(sym4_1, [32*32 length(sym4_1(1,1,:))]);
sym4_1 = sym4_1(ind_nnz,:);
subplot(2,4,2)
%mesh(mb42_1)
plot(t,sym4_1')
xline(examp,'--r');
xlabel('Time (sec)')
title('1D DWT SYM4')

%view(az,el)
sym4_1 = sym4_1 - min(sym4_1,[],'all');
sym4_1 = sym4_1 ./ (max(sym4_1,[],'all'));
e_overall_1 = entropy(sym4_1);
e_temporal_1 = 0;
for i = 1:548 %loop through each sensor
    e_temporal_1 = e_temporal_1 + entropy(sym4_1(i,:));
end
e_spatial_1 = 0;
for i = 1:49 %loop through each sensor
    e_spatial_1 = e_spatial_1 + entropy(sym4_1(:,i));
end

%(c) best 2d = db1
load('nmse_0.01_wavelet_db1_2D');
db1_2 = answers{2,1};
db1_2 = db1_2(:,:,startt:endt); %4900, 5030
examp_db1_2 = db1_2(:,:,exampi);
sq_db1_2 = zeros(224,224);
num = 1;
for i=1:7
    for j=1:7
        sq_db1_2(32*(i-1)+1:32*(i),32*(j-1)+1:32*(j)) = db1_2(:,:,num);
        num = num + 1;
    end
end
db1_2 = reshape(db1_2, [32*32 length(db1_2(1,1,:))]);
db1_2 = db1_2(ind_nnz,:);
subplot(2,4,3)
%mesh(db1_2)
plot(t,db1_2')
xline(examp,'--r');
xlabel('Time (sec)')
title('2D DWT DB1')

%view(az,el)
db1_2 = db1_2 - min(db1_2,[],'all');
db1_2 = db1_2 ./ (max(db1_2,[],'all'));
e_overall_2 = entropy(db1_2);
e_temporal_2 = 0;
for i = 1:548 %loop through each sensor
    e_temporal_2 = e_temporal_2 + entropy(db1_2(i,:));
end
e_spatial_2 = 0;
for i = 1:49 %loop through each sensor
    e_spatial_2 = e_spatial_2 + entropy(db1_2(:,i));
end

%(d) best 3d = mb4.2
load('nmse_0.01_wavelet_mb4.2_3D');
mb42_3 = answers{2,1};
mb42_3 = mb42_3(:,:,startt:endt); %4900, 5030
examp_mb42_3 = mb42_3(:,:,exampi);
sq_mb42_3 = zeros(224,224);
num = 1;
for i=1:7
    for j=1:7
        sq_mb42_3(32*(i-1)+1:32*(i),32*(j-1)+1:32*(j)) = mb42_3(:,:,num);
        num = num + 1;
    end
end
mb42_3 = reshape(mb42_3, [32*32 length(mb42_3(1,1,:))]);
mb42_3 = mb42_3(ind_nnz,:);
subplot(2,4,4)
plot(t,mb42_3')
xline(examp,'--r');
xlabel('Time (sec)')
title('3D DWT MB4.2')

mb42_3 = mb42_3 - min(mb42_3,[],'all');
mb42_3 = mb42_3 ./ (max(mb42_3,[],'all'));
e_overall_3 = entropy(mb42_3);
e_temporal_3 = 0;
for i = 1:548 %loop through each sensor
    e_temporal_3 = e_temporal_3 + entropy(mb42_3(i,:));
end
e_spatial_3 = 0;
for i = 1:49 %loop through each sensor
    e_spatial_3 = e_spatial_3 + entropy(mb42_3(:,i));
end

% plot the square images
subplot(2,4,5)
v = segment(:,:,exampi);
examp_original = reshape(v,[1024 1]);
examp_original = examp_original(ind_nnz);
v(logical(zero_zone)) = nan;
%surf(v)
imagesc(v)
xlabel('Taxel')
ylabel('Taxel')

subplot(2,4,6)
v = examp_sym4_1;
vv = reshape(v,[1024,1]);
vv = vv(ind_nnz);
err = nmse(examp_original,vv)
v(logical(zero_zone)) = nan;
%surf(v)
imagesc(v)
xlabel('Taxel')

subplot(2,4,7)
v = examp_db1_2;
vv = reshape(v,[1024,1]);
vv = vv(ind_nnz);
err = nmse(examp_original,vv)
v(logical(zero_zone)) = nan;
%surf(v)
imagesc(v)
xlabel('Taxel')

subplot(2,4,8)
v = examp_mb42_3;
vv = reshape(v,[1024,1]);
vv = vv(ind_nnz);
err = nmse(examp_original,vv)
v(logical(zero_zone)) = nan;
%surf(v)
imagesc(v)
xlabel('Taxel')

e_spatial = e_spatial / 548;
e_spatial_1 = e_spatial_1 / 548;
e_spatial_2 = e_spatial_2 / 548;
e_spatial_3 = e_spatial_3 / 548;

e_temporal = e_temporal / 49;
e_temporal_1 = e_temporal_1 / 49;
e_temporal_2 = e_temporal_2 / 49;
e_temporal_3 = e_temporal_3 / 49;

table = [e_temporal, e_spatial, e_overall; ...
    e_temporal_1, e_spatial_1, e_overall_1; ...
    e_temporal_2, e_spatial_2, e_overall_2; ...
    e_temporal_3, e_spatial_3, e_overall_3]';
%% Fig 10
load('filt_stag.mat')
load('ind_nnz.mat')
close all
startt = 1;
endt =  33550;




segment = filt_stag(:,:,startt:endt);
v_segment = reshape(segment, [32*32 length(segment(1,1,:))]);
vnnz_segment = v_segment(ind_nnz,:);

%plot(vnnz_segment')

load('nmse_0.01_wavelet_sym4_1D');
sym4_1 = answers{2,1};
sym4_1 = sym4_1(:,:,startt:endt); %4900, 5030
sym4_1 = reshape(sym4_1, [32*32 length(sym4_1(1,1,:))]);
sym4_1 = sym4_1(ind_nnz,:);

load('nmse_0.01_wavelet_db1_2D');
db1_2 = answers{2,1};
db1_2 = db1_2(:,:,startt:endt); %4900, 5030
db1_2 = reshape(db1_2, [32*32 length(db1_2(1,1,:))]);
db1_2 = db1_2(ind_nnz,:);

load('nmse_0.01_wavelet_mb4.2_3D');
mb42_3 = answers{2,1};
mb42_3 = mb42_3(:,:,startt:endt); %4900, 5030
mb42_3 = reshape(mb42_3, [32*32 length(mb42_3(1,1,:))]);
mb42_3 = mb42_3(ind_nnz,:);

close all

temp = .05;
spatial = .08;
t = [1:33550]/7;

time = length(vnnz_segment(1,:));
% 1D
figure()
subplot(2,3,1)
%plot temporally
errs = [];
for i=1:548
    errs =[errs nmse(vnnz_segment(i,:),sym4_1(i,:))];
end

plot(errs)

x = [0 0 temp temp];
y = [mean(errs)-std(errs)/2 mean(errs)+std(errs)/2 mean(errs)+std(errs)/2 mean(errs)-std(errs)/2];
hold on
fill(x,y,'red','FaceAlpha',0.3,'EdgeColor','red')
yline(mean(errs), 'red','LineWidth',2)
text(5,.047,['\mu = ',num2str(mean(errs),1)],'FontSize',7,'FontWeight','bold')
text(5,.042,['\sigma = ',num2str(std(errs),1)],'FontSize',7,'FontWeight','bold')

mean(errs)
std(errs)
title('Temporal Error')
ylim([0 temp])
xlim([1 548])
xlabel('Sensor Number')
ylabel('NMSE')

subplot(2,3,4)
errs = [];
for i=1:time
    errs =[errs nmse(vnnz_segment(:,i),sym4_1(:,i))];
end

plot(t,errs)

x = [0 0 max(t) max(t)];
y = [mean(errs)-std(errs)/2 mean(errs)+std(errs)/2 mean(errs)+std(errs)/2 mean(errs)-std(errs)/2];
hold on
fill(x,y,'red','FaceAlpha',0.3,'EdgeColor','red')
yline(mean(errs), 'red','LineWidth',2)
text(5,.077,['\mu = ',num2str(mean(errs),1)],'FontSize',7,'FontWeight','bold')
text(5,.070,['\sigma = ',num2str(std(errs),1)],'FontSize',7,'FontWeight','bold')

mean(errs)
std(errs)
yline(mean(errs), 'red','LineWidth',2)
title('Spatial Error')
ylim([0 spatial])
xlim([1 max(t)])
xlabel('Time (sec)')
ylabel('NMSE')

%2D
subplot(2,3,2)
%plot temporally
errs = [];
for i=1:548
    errs =[errs nmse(vnnz_segment(i,:),db1_2(i,:))];
end
plot(errs)

x = [0 0 temp temp];
y = [mean(errs)-std(errs)/2 mean(errs)+std(errs)/2 mean(errs)+std(errs)/2 mean(errs)-std(errs)/2];
hold on
fill(x,y,'green','FaceAlpha',0.3,'EdgeColor','green')
yline(mean(errs), 'green','LineWidth',2)
text(5,.047,['\mu = ',num2str(mean(errs),1)],'FontSize',7,'FontWeight','bold')
text(5,.042,['\sigma = ',num2str(std(errs),1)],'FontSize',7,'FontWeight','bold')

mean(errs)
std(errs)
title('Temporal Error')
ylim([0 temp])
xlim([1 548])
xlabel('Sensor Number')

subplot(2,3,5)
errs = [];
for i=1:time
    errs =[errs nmse(vnnz_segment(:,i),db1_2(:,i))];
end
plot(t,errs)

x = [0 0 max(t) max(t)];
y = [mean(errs)-std(errs)/2 mean(errs)+std(errs)/2 mean(errs)+std(errs)/2 mean(errs)-std(errs)/2];
hold on
fill(x,y,'green','FaceAlpha',0.3,'EdgeColor','green')
yline(mean(errs), 'green','LineWidth',2)
text(5,.077,['\mu = ',num2str(mean(errs),1)],'FontSize',7,'FontWeight','bold')
text(5,.070,['\sigma = ',num2str(std(errs),1)],'FontSize',7,'FontWeight','bold')

mean(errs)
std(errs)
title('Spatial Error')
ylim([0 spatial])
xlim([1 max(t)])
xlabel('Time (sec)')

%3D
subplot(2,3,3)
%plot temporally
errs = [];
for i=1:548
    errs =[errs nmse(vnnz_segment(i,:),mb42_3(i,:))];
end
plot(errs)

x = [0 0 temp temp];
y = [mean(errs)-std(errs)/2 mean(errs)+std(errs)/2 mean(errs)+std(errs)/2 mean(errs)-std(errs)/2];
hold on
fill(x,y,'blue','FaceAlpha',0.3,'EdgeColor','blue')
yline(mean(errs), 'blue','LineWidth',2)
text(5,.047,['\mu = ',num2str(mean(errs),1)],'FontSize',7,'FontWeight','bold')
text(5,.042,['\sigma = ',num2str(std(errs),1)],'FontSize',7,'FontWeight','bold')

mean(errs)
std(errs)
title('Temporal Error')
ylim([0 temp])
xlim([1 548])
xlabel('Sensor Number')

subplot(2,3,6)
errs = [];
for i=1:time
    errs =[errs nmse(vnnz_segment(:,i),mb42_3(:,i))];
end
plot(t,errs)

x = [0 0 max(t) max(t)];
y = [mean(errs)-std(errs)/2 mean(errs)+std(errs)/2 mean(errs)+std(errs)/2 mean(errs)-std(errs)/2];
hold on
fill(x,y,'blue','FaceAlpha',0.3,'EdgeColor','blue')
yline(mean(errs), 'blue','LineWidth',2)
text(5,.077,['\mu = ',num2str(mean(errs),1)],'FontSize',7,'FontWeight','bold')
text(5,.070,['\sigma = ',num2str(std(errs),1)],'FontSize',7,'FontWeight','bold')

mean(errs)
std(errs)
title('Spatial Error')
ylim([0 spatial])
xlim([1 max(t)])
xlabel('Time (sec)')
%% Fig 11
clear all
close all
load('filt_stag.mat')
load('ind_nnz.mat')
v_segment = reshape(filt_stag, [32*32 length(filt_stag)]);
vnnz_segment = v_segment(ind_nnz,:);% - 512;

all_peaks = [];
for i = 1:548
    signal = vnnz_segment(i,:);
    [pks,locs] = findpeaks(signal,'MinPeakProminence',20,'MinPeakDistance',50);
    
    % find all peaks
    peaks = [];
    for j = 1:length(locs)
        if locs(j)-25 < 1
            continue
        end
        if locs(j)+24 > 33550
            continue
        end
        peaks = [peaks; signal(locs(j) - 25 : locs(j) +24)];
    end
    all_peaks = [all_peaks; peaks];

end

plot(all_peaks')
figure()
% This is the average pressure! it looks like the best wavelets?
t = (1:50)/7;
plot(t,mean(all_peaks),'blue','LineWidth',10)
hold on
errorbar(t,mean(all_peaks),std(all_peaks),'.b')
%hold on
xlim([t(1) t(end)])
xlabel('Time (sec)')
ylabel('ADC Reading')
legend('Mean Tactile Interaction','Standard Deviation')

iter=5;
ysize = 300;
xsize = 100;

[phi,psi,xval] = wavefun('sym4',iter);
f=figure()
subplot(2,1,1)
plot(phi)
xlim([1 length(phi)])
title('Scaling Function')
subplot(2,1,2)
plot(psi)
xlim([1 length(psi)])
title('Wavelet Function')
sgtitle('Sym4')
f.Position = [100 100 xsize ysize];

[phi,psi,xval] = wavefun('bior4.4',iter);
f=figure()
subplot(2,1,1)
plot(phi)
xlim([1 length(phi)])
title('Scaling Function')
subplot(2,1,2)
plot(psi)
xlim([1 length(psi)])
title('Wavelet Function')
sgtitle('Bior4.4')
f.Position = [100 100 xsize ysize];

[phi,psi,xval] = wavefun('bior6.8',iter);
f=figure()
subplot(2,1,1)
plot(phi)
xlim([1 length(phi)])
title('Scaling Function')
subplot(2,1,2)
plot(psi)
xlim([1 length(psi)])
title('Wavelet Function')
sgtitle('Bior6.8')
f.Position = [100 100 xsize ysize];
%% Fig 12
clear all
close all
load('master_data_fixed.mat')
MSEs = [0.01,0.0043,0.0015,0.0001];

family_names = {'beyl','bior','bl','coif','db','dmey','fk','han','mb', ...
    'rbio','sym','vaid'};
shape = {"+","pentagram","v","square","o","x","^","<",">","hexagram", ...
    "diamond","*"}; % family name
N=15;
hex = ['#1f77b4';'#aec7e8';'#ff7f0e';'#ffbb78';'#2ca02c';'#98df8a';'#d62728';'#ff9896';'#9467bd';'#c5b0d5';'#8c564b';'#c49c94';'#e377c2';'#f7b6d2';'#7f7f7f';'#c7c7c7';'#bcbd22';'#dbdb8d';'#17becf';'#9edae5'];
raw = sscanf(hex.','#%2x%2x%2x',[3,Inf]).';
num = size(raw,1);
facecolors = raw(1+mod(0:N-1,num),:) / 255;
edgecolors = {'red','green','blue'}; %1D, 2D, 3D


f = figure()
Sparsities = [0;cell2mat(master_data(2:end,3))];


all_f = [];
all_s = [];
for j = 1:4
    subplot(2,2,5-j)
    fsize = [];
    spars = [];
    for i=2:900 %900
        % do only 1D data
        if (master_data{i, 10} == MSEs(j)) %if MSE = 0.01
            if (master_data{i, 9} == 1)
                % remove outlier rbio3.1 and dmey
                if (strcmp([master_data{i,6},num2str(master_data{i,7})], 'rbio3.1'))
                    continue
                end
                if (strcmp([master_data{i,6},num2str(master_data{i,7})],'dmey'))
                    continue
                end
                if (strcmp([master_data{i,6},num2str(master_data{i,7})],'coif5'))
                    continue
                end
                %if family name = ... set the shape (6)
                %if filter size = ... set the face color (8)
                %if dimension = ... set the edge color (9)
                x = strmatch(master_data{i, 6}, family_names);

                s = Sparsities(i);
                wname = [master_data{i,6},num2str(master_data{i,7})];
                [h,~,~,~] = wfilters(wname);
                l = length(h);

                fsize = [fsize, l];
                spars = [spars, s];
                plot(l, s, shape{x}, ...
                    'MarkerSize',17, 'MarkerEdgeColor', ...
                    edgecolors{master_data{i, 9}}, 'MarkerFaceColor', ...
                    facecolors(master_data{i, 8},:),'LineWidth', 1);
                % marker size = 20
                hold on
                %nump = nump+1;
            end
        end
    end
    % for each value in fsize, 
    [sortedF,I] = sort(fsize);
    sortedS = spars(I);
    % for each value in sortedF, get the average values (and STD) or S
    Fs = unique(sortedF);
    mSs = zeros(1,length(Fs));
    sdSs = zeros(1,length(Fs));
    for K = 1:length(Fs)
        ind = find(sortedF == Fs(K));
        mSs(K) = median(sortedS(ind));
        sdSS(K) = std(sortedS(ind));
    end


    %p = polyfit(Fs,mSs,2);
    [p,SD] = polyfit(sortedF,sortedS,1);
    x1 = linspace(1 , 24);
    [y1,delta] = polyval(p,x1,SD);
    plot(x1,y1,'LineWidth',10)
    plot(x1,y1+2*delta,'m--',x1,y1-2*delta,'m--')
    %legend('Linear Fit','95% Prediction Interval')

    all_f = [all_f; fsize];
    all_s = [all_s; spars];
    set(gca,'yscale','log')
    ylim([-inf max(spars)*1.1])
 %   set(gca,'xscale','log')
    xlabel('Filter Size')
    ylabel('Sparsity')
    title(['NMSE = ',num2str(MSEs(j))])
    f.Position = [100 100 300 400];
end