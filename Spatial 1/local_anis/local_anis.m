clc;clear;
n=256;
a=20; b=15;
field1= aniso_exp(n,a,b);
 imagesc(field1);colorbar
% figure();
% hist(reshape(field1,n^2,1),100)
% field1= aniso_exp_new(n,a,b);
% imagesc(field1);colorbar
fs=0.2; fm=3;
eps=7.5;
xcen=128; ycen=128;
c=10;
a=2;
H=fspecial('gaussian',4,2);
pic=iso_Exp_old(n,c,a);
img=stainmake(pic,field1,fm,fs,eps,xcen,ycen,H);
%  figure();
% hist(reshape(img,n^2,1),100)
%  H = fspecial('disk',3);
% H = fspecial('motion',5,0);
% 
% new=conv2(img,H);
% imagesc(new);colormap jet;colorbar
% figure();
% hist(reshape(new,numel(new),1),100)