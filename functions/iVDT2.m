function out = iVDT2(Data,ord1,ord2 )
sz=size(Data);
disp(sz);
lg=length(sz)-1;
% tmsz=sqrt(sz(1:end-1));
td=zeros(1,2*lg);
td(1:2:end)=ord1;
td(2:2:end)=ord2;

Data1=reshape(Data, [td sz(end)]);
% ld=zeros(2*lg);
% ld(1:2:end)=1:lg;
% ld(2:2:end)=lg+1:2*lg;
Data2=permute(Data1, [1:2:2*lg 2:2:2*lg 2*lg+1]);
sz2=size(Data2);
out= reshape(Data2, [prod(sz2(1:lg)) prod(sz2(lg+1:2*lg)) sz2(end)]);