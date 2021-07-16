clear;clc;
load ROI_Lin;
n = 256;
c=10; a=2;
G=d_maker(n,c,a);
Ga = fft2(G);
sim=100;
mu=zeros(64,sim); % # of ROI=64
for ii=1:sim
 xx1=iso_Exp(Ga,n);

  for j=1:64
  asdf=xx1(ROIs(:,:,j));
 mu(j,ii)=mean(asdf);
      end
end
mu0=mean(mu,2);

clearvars -except mu0;
save ('mu0','mu0');
