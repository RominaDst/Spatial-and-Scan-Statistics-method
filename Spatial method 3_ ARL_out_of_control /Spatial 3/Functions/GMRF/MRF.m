function [y1]= MRF(m,mucond,C, idx, D, condk,yk);
 y1=zeros(m,m);
 z=zeros(length(idx),1);
 x=zeros(length(idx),1);
% mucond=mu(idx)+(Qcond\(-D(idx,condk)*(yk-mu(condk))));
% C =chol(Qcond,'lower'); % Factor it for generation  %slow
% C=sparse(C);
z = randn(m^2-length(condk),1); % z is N(0,I)
x = mucond+C'\z;
y=zeros(m^2,1);
y(idx)=x;
y(condk)=yk;
y1=reshape(y,m,m);
% imagesc(y1)
% colorbar
end