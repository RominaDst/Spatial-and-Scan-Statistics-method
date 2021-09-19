clc;clear;
load beta_sem;load sig_sem; % load in-control historical data
sig1=log(sig);
s=mean(sig1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
b=mean(beta);
tmp=zeros(length(b),length(b));
S2=zeros(length(b),length(b)); % build covarince of in-control data
l=length(beta);
for k=1:l-1
    tmp=(1/(2*(l-1)))*(beta(k+1,:)-beta(k,:))'*(beta(k+1,:)-beta(k,:));
    S2=S2+tmp;
end
iS2=inv(S2); % covaince inverse
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n =32; % sumulated usrface size n*n
Z=zeros(n,n);
z=zeros(n^2,1);
[X,Y] = meshgrid(1:n,1:n); % construct gridded surface
yc=reshape(X,n^2,1);
xc=reshape(Y,n^2,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[~,W,~] = xy2cont(xc,yc); % construct a spatial contiguity matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x=[ones(n^2,1) xc yc]; %  autocorrelation test
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% general spatial model  
info.lflag  = 1;
beta=zeros(1,4); % estimated parameter betha
nos=1;% number of simulations
UCL=[22 ]; % upper control limit
xcen=n/2; ycen=n/2; % ceneter of simulated surface
fault=[1]; % numeber of fault
scenario=1; % shift scenario number
ARL=zeros(1,size(fault,2)); % average run lenght
xx=zeros(n,n);
k=3;
ucl=s+k*sqrt(var(sig1));% upper control limit for Shewhart chart of sigma
lcl=s-k*sqrt(var(sig1));% lower control limit for Shewhart chart of sigma
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for j=1:size(fault,2)
     counter=zeros(1,nos);
   for i=1:nos
       flag =0; 
       formatSpec = 'fault %3.0f and sim %3.0f \n';
       fprintf(formatSpec,j,i)
         while flag == 0
                 cd Functions
                 Z=all_maker(scenario); % simulated out of control for shift scenario 
                 cd ../
               z=reshape(Z,n^2,1);%image generating
              %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
              res = sem(z,x,W); % estimation of the spatial error model (SEM)
              beta=[res.beta(1) res.beta(2) res.beta(3) res.rho ];% construct beta for estimated parameters
              sig=res.sige; % model estimated error
               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                T1=(beta-b)*iS2* (beta-b)'; % Hotelling statistics for 
                T2=log(sig); % monitioring statistics for variance
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
               counter(i)=counter(i)+1;
                   if T1 > UCL || T2 > ucl || T2 < lcl
                    flag=1;
                    end
   end
          
  end
   
    ARL(1,j)=sum(counter)/nos;
    
end

ARL


