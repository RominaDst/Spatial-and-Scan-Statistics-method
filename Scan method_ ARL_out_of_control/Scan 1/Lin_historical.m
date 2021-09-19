% save the mean of  monitoring statistics of in-control historical data
clear;clc;
load ROI_Lin; % the locations of ROI in the gridded surface
n = 256; % size of the surface
c=10; a=2; % parameters of in-control surface
cd Functions
G=d_maker(n,c,a); % generate covarince matrix
Ga = fft2(G);
sim=10000; % number of simulation 
mu=zeros(64,sim); % # of ROI size = 64
for ii=1:sim
 xx1=iso_Exp(Ga,n);
  for j=1:64
  asdf=xx1(ROIs(:,:,j));
 mu(j,ii)=mean(asdf);
      end
end
mu0=mean(mu,2); mean of  monitoring statistics

clearvars -except mu0;
save ('mu0','mu0');
