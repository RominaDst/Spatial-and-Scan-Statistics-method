clc;clear;
n =32;
sim=100;
c=10; a=2;
G=d_maker(n,c,a);
Ga = fft2(G);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Z=iso_Exp(Ga,n);
z=reshape(Z,n^2,1);
[X,Y] = meshgrid(1:n,1:n);
yc=reshape(X,n^2,1);
xc=reshape(Y,n^2,1);
beta=zeros(sim,4);
sig=zeros(sim,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[~,W,~] = xy2cont(xc,yc);
W2 = slag(W,2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x=[ones(n^2,1) xc yc]; %  autocorrelation test
res1 = moran(z,x,W);
prt(res1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% general spatial model  

for i=1:sim
Z=iso_Exp(Ga,n);
z=reshape(Z,n^2,1);
   res = sem(z,x,W); 
%    prt(res)
%  zz=res.resid;
%  res1=sem(zz,x,W);
%  prt(res1);
%  
% res22 = moran(zz,x,W);
% prt(res22);
%  res =  sac(z,x,W,W2);
 beta(i,:)=[res.beta(1) res.beta(2) res.beta(3) res.rho ];
% beta(i,:)=[res.beta(1) res.beta(2) res.beta(3) res.rho res.lam];
%  beta(i,:)=[res.beta(1) res.beta(2) res.beta(3) res.rho res1.rho];
sig(i,1)=res.sige;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars -except beta sig;
save ('beta_sem','beta');
 save ('sig_sem','sig');
