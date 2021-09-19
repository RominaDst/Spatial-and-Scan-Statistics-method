%Gaussian stationary exponential spatial process
  function [cover]= wang_exp(r,c0,ce,ae,n)
ns=r^2;
% nm1=r-1;
% sig=1;
index=[1:(n/r):n];
x=index;
x=repmat(x,1,r);
y=index;
y=repmat(y,r,1);
y=reshape(y,1,numel(y));
cover=eye(ns,ns)*(ce);
counter=0;
for i=1:ns
    for j=i+1:ns
        counter=counter+1;
        h=sqrt((x(i)-x(j))^2+(y(i)-y(j))^2);
            cover(i,j)=c0+ce*(exp(-h/ae));
    end
end
cover=cover+cover';
cover=cover-ce*eye(ns);


