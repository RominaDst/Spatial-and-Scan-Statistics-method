clear;clc;
n = 256;% size of simulated surface
r=16; % sampling rate 
c0=0; % covarince parameters sill
theta=csvread('in_control_theta.csv');% in-control theta (sill and range of covarince)
theta1=mean(theta); % mean of in-control theta
xx=csvread('historical_data.csv'); % historical in-control surfaces
%%%%%%%%%%%%%%%
mu=mean(xx,2); %in-control mean
mur=reshape(mu,n,n);% reshape the mean
mu1=zeros(r,r);
index=[1:(n/r):n];
for i=1:r
  mu1(i,:)=mur(((i-1)*(n/r))+1,index);
end
mu1r=reshape(mu1,r^2,1); % reshape mean as a vector
%%%%%%%%%%%%%%%
cover=Grimshaw_exp(r,c0,theta1(1),theta1(2),n);% genarte the fitted covariance matrix
icover=inv(cover); % matrix inversion for Hoteeling statistics
%%%%%%%%%%%%%%%%%%%%%%%
UCL=319.4;% upper control limit
nos=1;% number os simulation
fault=[1]; % number of fault
ARL=zeros(1,size(fault,2)); % average run lenght
tmp1=zeros(r,r);
tmp2=zeros(r^2,1);
%%%%%%%%%%%%%%%%%%%%%%%
scenario=1; % shift scenatio
%%%%%%%%%%%%%%%%%%%%%%%
for j=1:size(fault,2)
    counter=zeros(1,nos);
   for i=1:nos
       flag =0;
        while flag == 0
              cd Functions
                 tmp=all_maker(scenario); %out_of_control surface under shift scenario
                 cd ../
                        for ii=1:r
                            tmp1(ii,:)=tmp(((ii-1)*(n/r))+1,index);
                        end
                   tmp2=reshape(tmp1,numel(tmp1),1);
                   T1=(tmp2-mu1r)' * icover * (tmp2-mu1r); % Hotelling statistics
                   counter(i)=counter(i)+1;
                   if T1 > UCL
                    flag=1;
                    end
        end
   end
   ARL(1,j)=sum(counter)/nos;
 
end
ARL

