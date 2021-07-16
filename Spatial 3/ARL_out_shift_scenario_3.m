clc;clear;
load beta_sem;load sig_sem;
sig1=log(sig);
s=mean(sig1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
b=mean(beta);
tmp=zeros(length(b),length(b));
S2=zeros(length(b),length(b));
l=length(beta);
for k=1:l-1
    tmp=(1/(2*(l-1)))*(beta(k+1,:)-beta(k,:))'*(beta(k+1,:)-beta(k,:));
    S2=S2+tmp;
end
iS2=inv(S2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n =32;
c=10; aa=2;
G=d_maker(n,c,aa);
Ga = fft2(G);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Z=zeros(n,n);
z=zeros(n^2,1);
[X,Y] = meshgrid(1:n,1:n);
yc=reshape(X,n^2,1);
xc=reshape(Y,n^2,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[~,W,~] = xy2cont(xc,yc);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x=[ones(n^2,1) xc yc]; %  autocorrelation test
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% general spatial model  
info.lflag  = 1;
beta=zeros(1,4);
nos=10000;
UCL=[22 ];
f=[0,0.25,0.5,1,1.5,2];
d=[4 2;-1 -0.5];
fault=combvec(f,d);
ARL=zeros(1,size(fault,2));
xx=zeros(n,n);
k=3;
ucl=s+k*sqrt(var(sig1));
lcl=s-k*sqrt(var(sig1));


tic
for j=1:size(fault,2)
     cd GMRF
     [idx,D,condk,mu,C,Qcond]= GMRF_input(n,fault(2,j),fault(3,j));
     yk=fault(1,j).*ones(n,1);
     mucond=mu(idx)+(Qcond\(-D(idx,condk)*(yk-mu(condk))));
    cd ../
     counter=zeros(1,nos);
   for i=1:nos
       flag =0; 
       formatSpec = 'fault %3.0f , sim %3.0f \n';
       fprintf(formatSpec,j,i)
         while flag == 0

                 cd GMRF
                 Z=MRF(n,mucond,C, idx, D, condk,yk);
                 cd ../
               z=reshape(Z,n^2,1);%image generating
              %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
              res = sem(z,x,W);
               beta=[res.beta(1) res.beta(2) res.beta(3) res.rho ];
                sig=res.sige;
               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                T1=(beta-b)*iS2* (beta-b)';
                T2=log(sig);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
               counter(i)=counter(i)+1;
                   if T1 > UCL || T2 > ucl || T2 < lcl
                
                    flag=1;
                    end
   end
          
  end
   

    ARL(1,j)=sum(counter)/nos;
    
end
toc


