clc;clear;
m = 256;
d1=4;
d2=-1;
GMRF_input(m,d1,d2);
load input
f=2; % mean conditional
yk=f.*ones(m,1);
mucond=mu(idx)+(Qcond\(-D(idx,condk)*(yk-mu(condk))));
tic
[y1]= MRF(m,mucond,C, idx, D, condk,yk);
toc
imagesc(y1);colorbar
