clear all;
close all;
addpath functions
addpath datasets
addpath tensorlab
addpath('Quality_Indices')

SNRm=25;
SNRh=25;

size_kernel=[8 8];
S=imread('original_rosis.tif');
F=load('R.mat');
S=double(S);


S=S(1:256,1:256,1:end);
% S=S/max(S(:));
F=F.R;
T=F(:,1:end);
 for band = 1:size(T,1)
        div = sum(T(band,:));
        for i = 1:size(T,2)
            T(band,i) = T(band,i)/div;
        end
 end
 
ratio=8;
s0=ratio/2;


%%%HSI
B=fspecial('average',size_kernel);
I_HSn = gauss3filter(S,B);
HSI=I_HSn(s0:ratio:end,s0:ratio:end,:);

%%%MSI
I_temp= reshape(S,size(S,1)*size(S,2),size(S,3));
I_ms=T*I_temp';
MSI=reshape(I_ms',size(S,1),size(S,2),size(T,1));

%%% add noise
if SNRm
    sigmam = sqrt(sum(MSI(:).^2)/(10^(SNRm/10))/numel(MSI));
    MSI = MSI + sigmam*randn(size(MSI));
end
if SNRh
    sigmah = sqrt(sum(HSI(:).^2)/(10^(SNRh/10))/numel(HSI));
    HSI = HSI + sigmah*randn(size(HSI));
end

%%%PaviaU
v1= [8 8 2 2];   
v2= [8 8 2 2];   
tv1 = [1 8 2 2];
tv2 = [1 8 2 2];
para.r=[0, 37, 3, 3, 2;
        0, 0,  9, 4, 2;
        0, 0,  0, 4, 2;
        0, 0,  0, 0, 2;
        0, 0,  0, 0, 0];

tic

para.max_iter=480;
para.mu=120;
para.lambda=0.1;
para.beta=0.1;
para.sigma=10;
para.v1=v1;
para.v2=v2;
[Lap,W]=Cp_Lap_spectral(HSI);
[Z,out,Utr,Ttr]=Fusion_FCTN(VDT2(MSI,v1,v2),VDT2(HSI,tv1,tv2),para,T,Lap,S);
quality_HCTR = QualityIndices(Z,S,ratio);
toc;
time_HCTR=toc;
