clc;clear;warning off all;fclose all;close all
%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\Parameters////////////////////////////////%
load('ROIs') % Locations of each ROI in the surdace , ROI size is 22*22, total number of ROIs =1580
load('Sim_Baseline') % contain the in-control parameters (mean and varince of each ROIs) 
pixs=[250;250]; % simulated surface size
m=10; %History Window
lratconst=(1./vars)/2; % the constant part of monitoring statistic R(m,s)
nums=X_count; % number of pixels in each ROIs
%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\///////////////////////////////////%
sim=21;
UCL=[7300];
fault=[1];
ARL=zeros(1,size(fault,2));
I=zeros(pixs(1),pixs(2));
xx1=zeros(pixs(1),pixs(2));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
scenario=2;

for j=1:size(fault,2)
    counter=zeros(1,sim-m);
    X_sum=zeros(1580,m);
                  for picc=1:m                      
                 cd Functions
                 I=all_maker(scenario);
                 cd ../
                    for p=1:1580
                          X_sum(p,picc)=sum(I(ROIs(:,:,p)));% sum of the pixels for each ROIs
                        end
                  end
     
    rat_tmp=zeros(1,sim-m);
    X_sum_new=zeros(1580,m);
    muhat=zeros(1580,1); % estimated mean 
    rats=zeros(1580,1);
    tt=zeros(1580,1);

for pic= m+1:sim    
    flag=0;
    while flag == 0
                 cd Functions
                 I=all_maker(scenario);
                 cd ../
        [tt]=evalROI(ROIs,I); % evaluate ROIs for calcualting the mean of 
        X_sum_new(:,1:m-1)=X_sum(:,2:m);
        X_sum_new(:,m)=tt(:,1);
        X_sum=X_sum_new;
        [rat_max]= maxim(m,X_sum,mus,nums,lratconst); %find the loaction of maximum R(s) in the past m images
        rat_tmp=(max(rat_max)); % monitioring statistics (R(m,s))
        counter(pic-m)=counter(pic-m)+1;
                   if rat_tmp > UCL
                    flag=1;
                   end
    
    end

  end
     ARL(1,j)=sum(counter)/(sim-m);
end

ARL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
