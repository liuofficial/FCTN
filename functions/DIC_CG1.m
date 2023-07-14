function [ E ] = DIC_CG1(E0,Z,O,Y,H,R,Lap,para)

R1=R'*R; 
O1=O*O';
H1=H*H';

H2=R'*Z*O'+para.lambda*Y*H';

n = length(H2(:));
maxiters = max(100,sqrt(n));
normb = norm(H2(:));

 r0=H2-R1*E0*O1-para.lambda.*E0*H1-para.beta.*Lap*E0;
 p0=r0;                       
 rtr=r0(:)'*r0(:);
 niters=0;
while sqrt(rtr)/normb > 2e-6  &&  niters < maxiters
    niters = niters+1;
    pp= R1*p0*O1+para.lambda.*p0*H1+para.beta.*Lap*p0;      
    pp1=p0(:)'*pp(:);         

    a=(r0(:)')*r0(:)/pp1;      
         E=E0+a*p0;            
         r1=r0-a*pp;           
         
         rtr=r1(:)'*r1(:);
         b1=(r1(:)'*r1(:))/(r0(:)'*r0(:)); 
         p1=r1+b1*p0;          
     p0=p1;                    
     r0=r1;                    
     E0=E;                    
 end
 
 
