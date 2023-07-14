function [Z,out,Utr,Ttr]=Fusion_FCTN(MSI,HSI,para,R,Lap,S)

max_psnr=0;
min_sam=0;
min_ergas=0;
max_q2n=0;
p_iter=0;
min_rmse=20000;
m_iter=0;

N=ndims(MSI);
Mway=size(MSI);
Hway=size(HSI);
Zway=Hway;
Zway(1)=Mway(1);
Mdim=para.r+para.r'+diag(Mway);
Hdim=para.r+para.r'+diag(Hway);
Sdim=Hdim;
Sdim(1,:)=Mdim(1,:);

Utr=cell(1,N);
Ttr=cell(1,N);
out=cell(1,N);
for i=1:N
    Utr{i}=rand(Mdim(i,:));
    Ttr{i}=rand(Hdim(i,:));
end

Ttr(2:end-1)=Utr(2:end-1);
Utr{end}=tmprod(Ttr{end},R,N);


for iter=1:para.max_iter 

    %%
    Y1=my_Unfold(HSI,Hway,1);
    tempH=tnreshape(tnprod_rest(Ttr,1),N,1);
    tempA=para.lambda.*tempH*tempH'+para.mu.*eye(size(tempH,1));
    Ttr{1}=my_Fold(para.lambda.*Y1*tempH'*pinv(tempA),Hdim(1,:),1);


    %%update U1
    Z1=my_Unfold(MSI,Mway,1);
    tempO=tnreshape(tnprod_rest(Utr,1),N,1);
    tempB=tempO*tempO'+para.mu.*eye(size(tempO,1));
    Utr{1}=my_Fold(Z1*tempO'*pinv(tempB),Mdim(1,:),1);

    %%update Ut
    for t=2:N-1
        Zt=my_Unfold(MSI,Mway,t);
        Yt=my_Unfold(HSI,Hway,t);
        tempOt=tnreshape(tnprod_rest(Utr,t),N,t);
        tempHt=tnreshape(tnprod_rest(Ttr,t),N,t);
        tempA=Zt*tempOt'+para.lambda.*Yt*tempHt';
        tempB=tempOt*tempOt'+para.lambda.*tempHt*tempHt'+para.mu.*eye(size(tempHt,1));
        Utr{t}=my_Fold(tempA*pinv(tempB),Mdim(t,:),t);
        Ttr{t}=Utr{t};
    end

    %%Update Ud+1
    Zd=my_Unfold(MSI,Mway,N);
    Yd=my_Unfold(HSI,Hway,N);
    Od=tnreshape(tnprod_rest(Utr,N),N,N);
    Hd=tnreshape(tnprod_rest(Ttr,N),N,N);
    E0=my_Unfold(Ttr{N},Hdim(N,:),N);
    Ttr{end}=my_Fold(DIC_CG1(E0,Zd,Od,Yd,Hd,R,Lap,para),Hdim(N,:),N);
    Utr{end}=tmprod(Ttr{end},R,N);


    %%%
    out=Utr;
    out(end)=Ttr(end);
    out1=my_Unfold(out{1},Sdim(1,:),1);
    out2=tnreshape(tnprod_rest(out,1),N,1);
    Z=my_Fold(out1*out2,Zway,1);
    Z=iVDT2(Z,para.v1,para.v2);
    
    rmse=getrmse(double((S)),double((Z)));
    disp(['iter:',num2str(iter),',rmse:',num2str(rmse)]);
    if(rmse<min_rmse)
        min_rmse=rmse;
        m_iter=iter;
    end
    if(mod(iter,10)==1)
        quality_HCTR = QualityIndices(Z,S,8);
        if(quality_HCTR.psnr>=max_psnr)
            max_psnr=quality_HCTR.psnr;
            min_sam=quality_HCTR.sam;
            min_ergas=quality_HCTR.ergas;
            p_iter=iter;
        end
         disp(['max_psnr:',num2str(max_psnr),',min_sam:',num2str(min_sam),',min_ergas:',num2str(min_ergas)']);
    end
end
disp(['min_rmse:',num2str(min_rmse),',m_iter:',num2str(m_iter),'max_psnr:',num2str(max_psnr),'min_sam:',num2str(min_sam),',min_ergas:',num2str(min_ergas),',p_iter:',num2str(p_iter)]);



end


