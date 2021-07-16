clc;clear;close all;
t1=csvread('wang_theta_historical.csv');
t1=t1(:,2:3);
t=mean(t1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tmp=zeros(2,2);
S2=zeros(2,2);
l=length(t1);
for k=1:l-1
    tmp=(1/(2*(l-1)))*(t1(k+1,:)-t1(k,:))'*(t1(k+1,:)-t1(k,:));
    S2=S2+tmp;
end
iS2=inv(S2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n = 2^8;c=10;a=2;
G=d_maker(n,c,a);
Ga = fft2(G);
Z= iso_Exp(Ga,n);
[resy0,resx0]=size(Z);
[p0,n0]=variogram_maker1(resx0);
np=n0(2:resx0/2+1);
hs=[1:resx0/2];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
r=16;
index=[1:(n/r):n];
x=index;
x=repmat(x,1,r);
y=index;
y=repmat(y,r,1);
y=reshape(y,1,numel(y));
B=[ones(r^2,1) x' y'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% options = optimoptions(@fmincon,'GradObj','on','Display','iter');
options = optimoptions(@fmincon,'GradObj','on');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 theta=zeros(1,3);
 flag=0;
% counter2=0;
 while flag==0
%      counter2=counter2+1;
                clf
                nos=500;
                for i=1:nos
               Z= iso_Exp(Ga,n);  %image generating
%               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
              [gam]=emprical_variogram(Z,resy0,resx0,p0,n0);  %estimating parameters
              fun=@(x)fun2(x,np,gam,hs);
              theta= fmincon(fun,[.5,1,12],[],[],[],[],[0 0 0],[],[],options);
             for j=1:  length(hs)
                  var(j)= theta(1)+ theta(2)*(1-exp(-hs(j)/theta(3)));
              end
              hold on
              plot(hs,var,'r')
              
              T2(i)=(theta(2:3)-t)*iS2*(theta(2:3)-t)';
%               ttmp(counter2)=T2;
                end
                
                
                
              if T2>20
%                   plot(hs,var,'r')
                  break
%               else
%                   plot(hs,var,'k')
%                   pause(1)
%                   clc
%               end
              end       
 pause(3)
 end
 
plot(T2,'*');grid on
m=200;
p=2;
ucl1=finv(0.995,p,m-p)*(p*(m+1)*(m-1))/(m*(m-p))
ans=sum(T2>ucl1)/nos
ans=sum(T2>15.5)/nos
1/ans


%     ARL(j)=sum(counter)/nos;
%     MRL(j)=median(counter);
%   end
%  
%   save( 'ARL2_incontrol.mat', 'ARL','UCL','nos','MRL'); 
 
   plot(T2,'*');grid on