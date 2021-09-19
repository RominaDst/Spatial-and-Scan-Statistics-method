% This code generates a two dimensional Gaussian process with mean zero and
% exponentional covarince with using circulant embedding method
clc;clear;
n = 256; % size of simulated surface n*n
c=10; a=2; % parameters of covarince function
G=covaraince_maker(n,c,a); % generating covaraince matrix
Ga = fft2(G);
Z = randn(n) + sqrt(-1)*randn(n);
X = real(fft2((Ga.^.5).*Z/n));
imagesc(X); colormap(bone);colorbar;

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

