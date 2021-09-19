clear;clc;
load ROI; % load the output of "ROI_maker.m"  : generated region of interest location of the 256*256 grid with defined scan_window
load mu0; % mean of in-control historical data : load the output of "Lin_historical.m"
n = 256;% size of simulated surface n*n
l=0.2; %spatial parameters
w=0.1; %temporal parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%scan window size is 32 * 32
sl=0; %spatially weighted constant
for i=1:16
    sl=sl+(1-l)^i;
end
xx1=zeros(n,n);
cen=[16 16]; % location the center of scan window
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
UCL=[ 3.345 ]; % upper control limit (number otained by simulation to set ARL_in_control =200)
nos=1; % number of simulation 
ARL=zeros(1,size(1,2)); % Average Run Lenght
tmp=zeros(n,n); % pre_allocated variabe for out_of_control_surface (saving memory)
xx=zeros(n,n); % pre_allocated variabe for in_control_surface (saving memory)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fault=[1]; % numebr of faults
scenario=3; % shift scenario numebr according to the paper {1,2,3,4,5}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ARL_simulation

for jj=1:size(fault,2)
    counter=zeros(1,nos); % Variable Run Length
   for i=1:nos
       flag =0;
       xx1=zeros(n,n);
        while flag ==0        
           %%%%%%%% Simulating out_of_control surface %%%%%%
           cd Functions           
           tmp=all_maker(scenario);
           cd ../
           % Use the selected shift scenario and copy the corresponding
           % section and paste here.
           % The variable "tmp" should be generated here and is the out_of_control surface for selected shift scenario.
           %%%%%%%%%%%%%%%%%%      
            xx1=(1-w).*xx1 + w.*(tmp); % EWMA statistic : E_((m,n) ) (t)=(1-w) E_((m,n) ) (t-1)+w(X_((m,n) ) (t)-mu_0
              S=zeros(1,64);
              k=0;
                 for j=1:64
                    ss=0;
                     asdf=xx1(ROIs(:,:,j)); % load region of interest 
                      asdf1=reshape(asdf,32,32); % reshape the ROI as a square
                              for d=1:16 
                                x=asdf1(cen(1)-d+1:cen(1)+d,cen(2)-d+1:cen(2)+d);
                                asdf1(cen(1)-d+1:cen(1)+d,cen(1)-d+1:cen(1)+d)=0; % calculate Chebyshev distance 
                                mu1=mean(mean(nonzeros(x)));
                                ss=(((1-l)^d)/sl)*(mu1-(mu0(j,:)))+ss; % weighted moving window
                              end
                              S(j)=ss;
                           end  
                         k=max(S)/sqrt(var(S)); % Monitoring statistics 
               counter(i)=counter(i)+1;
                   if k > UCL  
                    flag=1;
                   end
        end
   end
   ARL(jj)=sum(counter)/nos;
end

ARL