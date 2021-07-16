clc;clear;
b1=csvread('wang_beta_historical.csv');
t1=csvread('wang_theta_historical.csv');
t1=t1(:,2:3);
t=mean(t1);
b=mean(b1);
t=mean(t1);
%%%%%%%%%%%%%%
n = 2^8;c=10;a=2;
G=d_maker(n,c,a);
Ga = fft2(G);
Z= iso_Exp(Ga,n);
[resy0,resx0]=size(Z);
[p0,n0]=variogram_maker1(resx0);
np=n0(2:resx0/2+1);
hs=[1:resx0/2];
%%%%%%%%%%%%%%%
tmp=zeros(2,2);
S2=zeros(2,2);
l=length(t1);
for k=1:l-1
    tmp=(1/(2*(l-1)))*(t1(k+1,:)-t1(k,:))'*(t1(k+1,:)-t1(k,:));
    S2=S2+tmp;
end
% S2=cov(t1);
iS2=inv(S2);
%%%%%%%%%%%%%%
r=16;
index=[1:(n/r):n];
x=index;
x=repmat(x,1,r);
y=index;
y=repmat(y,r,1);
y=reshape(y,1,numel(y));
B=[ones(r^2,1) x' y'];
%%%%%%%%%%%%%%%%%%%%%
options = optimoptions(@fmincon,'GradObj','on');
nos=500;
% ARL=zeros(1,nos);
% MRL=zeros(1,nos);
% UCL=[8.9];
xcen=128;
ycen=128;
eps=2;
H=fspecial('gaussian',4,2);
%%%%%%%%%%%%%%%%%%%%%
% 
%  for j=1:length(UCL)
%     counter=zeros(1,nos);
   for i=1:nos
%        flag =0;
         fault=[2,0.05];
%          while flag == 0
              z1=zeros(r,r);
              theta=zeros(1,3);
              beta=zeros(1,3);
%               Z= iso_Exp(Ga,n);  %image generating
              cd stain
                 xx=iso_Exp(Ga,n);
                 Z=stainmake(xx,fault(1),fault(2),eps,xcen,ycen,H);
                    imagesc(Z)
               cd ../
              %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
              for jj=1:r
                z1(jj,:)=Z(((jj-1)*(n/r))+1,index); %downsampling
              end
              z=reshape(z1,r^2,1);
              %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
              [gam]=emprical_variogram_mex(Z,resy0,resx0,p0,n0);  %estimating parameters
              fun=@(x)fun2(x,np,gam,hs);
              theta= fmincon(fun,[.5,1,12],[],[],[],[],[0 0 0],[],[],options);
              sf = fit([x', y'],z,'poly11');
              beta= coeffvalues(sf);
              %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
              covert=wang_exp(r,theta(1),theta(2),theta(3),n);
              icovert=inv(covert);
              fisher1=B'*icovert*B;
               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
              T1=(beta-b)*fisher1*(beta-b)'
               T2=(theta(2:3)-t)*iS2*(theta(2:3)-t)'
               beta
               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                counter(i)=counter(i)+1;
%                    if T1 > UCL(j) 
%                     flag=1;
%                     end
%         end
  end
%    plot(T1,'*')
%    ARL(j)=sum(counter)/nos;
%    MRL(j)=median(counter);
% end
 
%  save( 'ARL1_incontrol.mat', 'ARL','UCL','nos','MRL'); 
plot(T1,'*');grid on
% m=200;
% p=3;
% ucl1=finv(0.995,p,m-p)*(p*(m+1)*(m-1))/(m*(m-p))
% (sum(T1>ucl1)/nos)
(sum(T2>22)/nos)

% ucl2=betainv(0.95,p/2,(m-p-1)/2)*((m-1)^2)/m


    