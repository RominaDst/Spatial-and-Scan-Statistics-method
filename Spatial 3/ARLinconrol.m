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
c=10; a=2;
G=d_maker(n,c,a);
Ga = fft2(G);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Z=zeros(n,n);
z=zeros(n^2,1);
[X,Y] = meshgrid(1:n,1:n);
yc=reshape(X,n^2,1);
xc=reshape(Y,n^2,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[~,W,~] = xy2cont(xc,yc);
W2 = slag(W,2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x=[ones(n^2,1) xc yc]; %  autocorrelation test
% res1 = moran(z,x,W);
% prt(res1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% general spatial model  
sim=1;
T1=zeros(1,sim);
T2=zeros(1,sim);
info.lflag  = 1;
eps=2;
xcen=16; ycen=16;
fault=[4;0.5];
H=fspecial('gaussian',4,2);
xx=zeros(n,n);

tic
for i=1:sim
  
%                cd stain
%                  xx=iso_Exp(Ga,n);
%                  fault=[2;0.1];
%                  Z=stainmake(xx,fault(1,1),fault(2,1),eps,xcen,ycen);
%                    imagesc(Z);colorbar
%                cd ../
%  Z=iso_Exp(Ga,n);
z=reshape(Z,n^2,1);
 res = sem(z,x,W);
  beta=[res.beta(1) res.beta(2) res.beta(3) res.rho ];
% beta=[res.beta(1) res.beta(2) res.beta(3) res.rho res.lam];
   sig=res.sige;
   T1=(beta-b)*iS2* (beta-b)'
  
   T2=log(sig);
   T2-ucl
end
toc
imagesc(Z);colorbar
% loc=find(T1>17);
% loc(1)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 plot(T1,'*');grid on
% m=100;
% p=length(b);
% ucl1=finv(0.995,p,m-p)*(p*(m+1)*(m-1))/(m*(m-p))
% sum(T1>ucl1)/sim
 1/(sum(T1>19.6)/sim)

% figure()
% plot(T2,'*');grid on
 k=3;
 ucl=s+k*sqrt(var(sig1))
 lcl=s-k*sqrt(var(sig1))
% (sum(t>=ucl)+sum(t<=lcl))/960
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
