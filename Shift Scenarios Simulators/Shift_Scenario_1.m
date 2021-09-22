% This code generates local mean shift.
%%%%%%%%%%%%%%%%
%in-control surface
clc;clear;
n = 256; % size of simulated surface n*n
c=10; a=2; % parameters of covarince function
G=covaraince_maker(n,c,a); % generating covaraince matrix
Ga = fft2(G);
Z = randn(n) + sqrt(-1)*randn(n);
X = real(fft2((Ga.^.5).*Z/n));
%imagesc(X); colormap(bone);colorbar;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%induced local mean shift
xcen=120; ycen=120; % location of the center of shift on the x and y axis
fm=[4]; % magnitude of the mean shift
fs=[0.09]; % shift size ratio (comapre to the original suface size n*n)
cd stain % The folder "stain" contains the NURBS codes for generating the random blot shaped shift
H=fspecial('gaussian',4,2); % create two-dimensional Guassian filter and will be applied on the border of the shift area
image_scenario_1=stainmake(X,fm,fs,eps,xcen,ycen,H); % function that induces the genearted shift on the in-control surface
imagesc(image_scenario_1); colormap(bone);colorbar;
%%%%%%%%%%%%%%%%%%%%Functions%%%%%%%%%%%%%%%%%%%%%%%%%

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

