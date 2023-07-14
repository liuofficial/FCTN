function [Lap W]=Cp_Lap_spectral(HSI)
% sz_hsi=sz(HSI);
% sz_msi=sz(MSI);

H_mat4=tens2mat(HSI,3);   
% M_mat4=tens2mat(MSI,4);
alpha=10;
% construct the Lapacian matrix
options = [];
      options.NeighborMode = 'KNN';
      options.k =20;
      options.WeightMode ='HeatKernel';% 'Binary';
      options.t =10;
      options.bSelfConnected=0;
      options2=options;
      options2.t=0.5;
%       W1= constructW(H_mat4,options);
%       W2= constructW(M_mat4,options2);
%       W=W1+W2;
W= constructW([H_mat4],options);
%  W=W./repmat(sum(W,2),1,size(W,2));
%       W=W./repmat(sum(W,2),1,size(W,2));
      Lap=diag(sum(W,2))-W;
