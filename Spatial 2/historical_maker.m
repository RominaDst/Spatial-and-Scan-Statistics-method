clc;clear;close all;
sim=200;
n = 2^8;c=10;a=2;
G=d_maker(n,c,a);
Ga = fft2(G);
theta=zeros(sim,3);
beta=zeros(sim,3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Z= iso_Exp(Ga,n);
[resy0,resx0]=size(Z);
[p0,n0]=variogram_maker1(resx0);
np=n0(2:resx0/2+1);
hs=[1:resx0/2];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[X,Y] = meshgrid(1:n,1:n);
x=reshape(X,n^2,1);
y=reshape(Y,n^2,1);
z=zeros(n^2,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
options = optimoptions(@fmincon,'GradObj','on');

for i=1:sim

Z= iso_Exp(Ga,n);
% [gam]=emprical_variogram_mex(Z,resy0,resx0,p0,n0);
% % gam=gam(1:64);
% % np=np(1:64);
% % hs=hs(1:64);
% fun=@(x)fun2(x,np,gam,hs);
% 
%  theta(i,:) = fmincon(fun,[.5,1,12],[],[],[],[],[0 0 0],[],[],options);
 
  z=reshape(Z,n^2,1);
  sf = fit([x, y],z,'poly11');
 
  beta(i,:) = coeffvalues(sf);
end

%  csvwrite('wang_theta_historical.csv',theta);
 csvwrite('wang_beta_historical.csv',beta);

