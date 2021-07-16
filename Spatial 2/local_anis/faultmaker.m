% close all;clc;clear
n = 256;
r=16; c0=0; c=10; a=2;
G=d_maker(n,c,a);
Ga = fft2(G);
xcen=220; ycen=220;
% fm=[1,2,3,4];
% fs=[0.05,0.06,0.07,0.08];
% % fault=combvec(fm,fs);
%  fault=[2,3,4,5,2,3,4,5,2,3,4,5,2,3,4,5;0.0400000000000000,0.0400000000000000,0.0400000000000000,0.0400000000000000,0.0500000000000000,0.0500000000000000,0.0500000000000000,0.0500000000000000,0.0700000000000000,0.0700000000000000,0.0700000000000000,0.0700000000000000,0.0800000000000000,0.0800000000000000,0.0800000000000000,0.0800000000000000];
% eps=2;
 fault=[3;0.06];
eps=2;
xx1=iso_Exp(Ga,n);
std(reshape(xx1,n^2,1))
% fault=[4;0.06]
%  H = fspecial('disk',3);
% H = fspecial('motion',5,0);
% fm=3;
% fs=0.06;
H=fspecial('gaussian',4,2);
for j=1: size(fault,2)
%  tmp=stainmake(xx1,fm,fs,eps,xcen,ycen,H);
 tmp=stainmake(xx1,fault(1,j),fault(2,j),eps,xcen,ycen,H);

% new=conv2(tmp,H);
 figure();imagesc(tmp);colorbar;colormap bone

end