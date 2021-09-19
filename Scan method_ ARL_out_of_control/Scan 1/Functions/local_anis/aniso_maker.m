function [field1]= aniso_maker(n,BlkCirc_row);
m=n;
lam=real(fft2(BlkCirc_row))/(2*m-1)/(2*n-1);
lam(lam(:)<0)=0; lam=sqrt(lam);
% generate field with covariance given by block circulant matrix
F=fft2(lam.*complex(randn(2*m-1,2*n-1),randn(2*m-1,2*n-1)));
F=F(1:m,1:n); % extract sub-block with desired covariance
field1=real(F);
end