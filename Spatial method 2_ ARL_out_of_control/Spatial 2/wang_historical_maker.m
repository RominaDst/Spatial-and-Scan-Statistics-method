clc;clear;close all;
sim=20;
t1=csvread('wang_theta_historical.csv');

t=mean(t1);
sig=cov(t1);
% options = optimoptions(@fmincon,'TolCon',0.01,'TolFun',0.01,'MaxFunEvals',10000,'MaxIter',10000);

%%%%%%%%%%%%%%
tmp=zeros(3,3);
S2=zeros(3,3);
for k=1:sim-1
    tmp=(1/(2*(sim-1)))*(t1(k+1,:)-t1(k,:))'*(t1(k+1,:)-t1(k,:));
%     v(k,:)=(t1(k+1,:)-t1(k,:))';
    S2=S2+tmp;
end
n = 2^8;c=10;a=2;
G=d_maker(n,c,a);
Ga = fft2(G);
theta=zeros(sim,3);
beta=zeros(sim,3);

Z= iso_Exp(Ga,n);
[resy0,resx0]=size(Z);
[p0,n0]=variogram_maker1(resx0);
np=n0(2:resx0/2+1);
hs=[1:resx0/2];

[X,Y] = meshgrid(1:n,1:n);
x=reshape(X,n^2,1);
y=reshape(Y,n^2,1);
T=zeros(sim,1);
iS2=inv(S2);
 tic
 sim=50;
for i=1:sim

Z= iso_Exp(Ga,n);

% z(:,i)=reshape(Z,numel(Z),1);
% imagesc(Z);colorbar
% tic
[gam]=emprical_variogram(Z,resy0,resx0,p0,n0); 
%  [gam]=emprical_variogram_mex(Z,resy0,resx0,p0,n0);
% toc

fun=@(x)fun2(x,np,gam,hs);
% fun=@(theta)fitvaro(theta,np,gam,hs);


% options = optimoptions(@fmincon,'Algorithm','active-set','ConstraintTolerance',0.01,'OptimalityTolerance',0.01);

%  theta = fmincon(fun,[0,1,12],[],[],[],[],[0 0 0],[]);
 
%  options = optimoptions(@fmincon,'SpecifyObjectiveGradient',true);
 options = optimoptions(@fmincon,'GradObj','on');
%  options = optimoptions(@fmincon,'Display','Iter','DiffMaxChange',.0001);
 [theta,FVAL,EXITFLAG,OUTPUT,LAMBDA,GRAD,HESSIAN] = fmincon(fun,[.5,1,12],[],[],[],[],[0 0 0],[],[],options)
 asdf(i,:)=theta;
%   toc
% % [x exitflag]=fmincon(fun,[0,1,12],[],[],[],[],[0 0 0],[])
   T(i,1)=(theta-t)*iS2*(theta-t)';
% if isnan(T(i,1))
%     asdf=3;
% end
%  sf = fit([x, y],z(:,i),'poly11');
% 
%  beta(i,:) = coeffvalues(sf);
end
 toc
 plot(T,'*')
m=50;
p=3;
ucl1=finv(0.95,p,m-p)*(p*(m+1)*(m-1))/(m*(m-p))
ucl2=betainv(0.95,p/2,(m-p-1)/2)*((m-1)^2)/m
%  csvwrite('wang_theta_historical.csv',theta);

ucl11=finv(0.5,p,m-p)*(p*(m-1)*(m-1))/(m*(m-p))
ucl21=betainv(0.5,p/2,(m-p-1)/2)*((m-1)^2)/m
