clc;clear;
b1=csvread('beta_historical.csv'); % load historical in-control beta
t1=csvread('theta_historical.csv'); % laod historical in-control theta
t1=t1(:,2:3);
t=mean(t1); % mean of in -control theta
b=mean(b1); % mean of in-control beta
%%%%%%%%%%%%%%
n = 2^8;% simulated surface size n*n
c=10;a=2; % in-control sufrace parameter
cd Functions
G=d_maker(n,c,a);
Ga = fft2(G);
Z= iso_Exp(Ga,n); % in-control surface
cd ..
[resy0,resx0]=size(Z);
[p0,n0]=variogram_maker1(resx0);% variogram of in-control surface
np=n0(2:resx0/2+1);
hs=[1:resx0/2];
%%%%%%%%%%%%%%%
tmp=zeros(2,2);
S2=zeros(2,2); % in-control covarince
l=length(t1);
for k=1:l-1
    tmp=(1/(2*(l-1)))*(t1(k+1,:)-t1(k,:))'*(t1(k+1,:)-t1(k,:));
    S2=S2+tmp;
end
iS2=inv(S2); % inverse of in-control covarince
%%%%%%%%%%%%%%
r=16; % sampling rate
index=[1:(n/r):n];
x=index;
x=repmat(x,1,r);
y=index;
y=repmat(y,r,1);
y=reshape(y,1,numel(y));
B=[ones(r^2,1) x' y']; % varaible B from joint distribution
%%%%%%%%%%%%%%%%%%%%%
options = optimoptions(@fmincon,'GradObj','on'); % optimization setting
nos=1; % simulation number
fault=[1]; % number of fault
UCL1=[17]; % upper contol limit for first monitoring statistics 
UCL2=[22]; % % upper contol limit for second monitoring statistics 
xx=zeros(n,n);
scenario=1; % shift scenario
%%%%%%%%%%%%%%%%%%%%%

 for j=1:size(fault,2)
    counter=zeros(1,nos);
   for i=1:nos
       flag =0;
         while flag == 0
              z1=zeros(r,r);
              theta=zeros(1,3);
              beta=zeros(1,3);
              cd Functions
              Z=all_maker(scenario); % generating out-of_control surface for shift scenario
               cd ../
              %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
              for jj=1:r
                z1(jj,:)=Z(((jj-1)*(n/r))+1,index); %downsampling
              end
              z=reshape(z1,r^2,1);
              %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
              [gam]=emprical_variogram(Z,resy0,resx0,p0,n0);  %estimating parameters from emorical variogram
              fun=@(x)fun2(x,np,gam,hs);
              theta= fmincon(fun,[.5,1,12],[],[],[],[],[0 0 0],[],[],options);% minimization function
              sf = fit([x', y'],z,'poly11'); % fit polynomial on mean surface
              beta= coeffvalues(sf); % estimated coefficients
              %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
              covert=wang_exp(r,theta(1),theta(2),theta(3),n); % Gaussian stationary exponential spatial process
              icovert=inv(covert);
              fisher1=B'*icovert*B; % Fisher information matrix 
               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
               T1=(beta-b)*fisher1*(beta-b)'; % Hotelling statistics for first varaibles
               T2=(theta(2:3)-t)*iS2*(theta(2:3)-t)'; % Hotelling statistics for second varaibles
               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
               counter(i)=counter(i)+1;
                    if T1 > UCL1 || T2> UCL2
                    flag=1;
                   end
        end
 end

   ARL(1,j)=sum(counter)/nos;
   
 end

ARL



    
