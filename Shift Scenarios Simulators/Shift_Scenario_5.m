% This code generates generates Local Anisotropic Stationary Correlation shift
%%%%%%%%%%%%%%%%
clc;clear;
n = 256; % size of simulated surface n*n
c0=0; c=10; aa=2; % parameters of in-control covarince function
a=20; b=15; %  parameters of Anisotropic  covarince function
G=covaraince_maker(n,c,aa); % generating covaraince matrix
Ga = fft2(G);
cd local_anis % folder contain codes for generating blot shaped shift with NURBS 
[BlkCirc_row]=rho_maker(n,a,b); %generate block circulant matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%induced local mean shift
xcen=120; ycen=120; % location of shift on the x and y axis
fm=[4]; % magnitude of the mean shift
fs=[0.2]; % shift size ratio (comapre to the original suface size n*n)
eps=7.5; % parameter to control the shape of the stain
H=fspecial('gaussian',4,2);% Guassian 2D Filter
field1= aniso_maker(n,BlkCirc_row); % % generate field with covariance given by block circulant matrix
pic=iso_Exp(Ga,n); % generate in-control surface
x=stainmake(pic,field1,fm,fs,eps,xcen,ycen,H); % combine loacl shift with anistropic shift
imagesc(x);colormap(bone)
%%%%%%%%%%%%%%%%%%%%Functions%%%%%%%%%%%%%%%%%%%%%%%%%

function [BlkCirc_row]=rho_maker(n,a,b);
m=n; 
% size of covariance matrix is m^2*n^2
tx=[0:n-1]; ty=[0:m-1]; % create grid for field

rho=@(x,y)((1-x^2/a^2-x*y/(b*a)-y^2/b^2)...
*exp(-(x^2/a^2+y^2/b^2)));

Rows=zeros(m,n); Cols=Rows;
for i=1:n
for j=1:m
Rows(j,i)=rho(tx(i)-tx(1),ty(j)-ty(1)); % rows of blocks
Cols(j,i)=rho(tx(1)-tx(i),ty(j)-ty(1)); % columns
end
end
% create the first row of the block circulant matrix
% with circulant blocks and store it as a matrix suitable for fft2
BlkCirc_row=[Rows, Cols(:,end:-1:2);
Cols(end:-1:2,:), Rows(end:-1:2,end:-1:2)];
end

function [field1]= aniso_maker(n,BlkCirc_row);
m=n;
lam=real(fft2(BlkCirc_row))/(2*m-1)/(2*n-1);
lam(lam(:)<0)=0; lam=sqrt(lam);
% generate field with covariance given by block circulant matrix
F=fft2(lam.*complex(randn(2*m-1,2*n-1),randn(2*m-1,2*n-1)));
F=F(1:m,1:n); % extract sub-block with desired covariance
field1=real(F);
end

function G=covaraince_maker(n,c,a);
t1 = [0:1/n:1-1/n];
t2 = t1;
for i=1:n 
    for j=1:n
        G(i,j)=exp(-c*(sqrt(min(abs(t1(1)-t1(i)), ...
        1-abs(t1(1)-t1(i)))^2 + min(abs(t2(1)-t2(j)), ...
        1-abs(t2(1)-t2(j)))^2)))^a;
    end;
  end;
end

