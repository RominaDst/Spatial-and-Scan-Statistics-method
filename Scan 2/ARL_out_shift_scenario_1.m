clc;clear;warning off all;fclose all;close all
%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\Parameters////////////////////////////////%
load('ROIs')
load('Sim_Baseline') % historical in-control data that needs to be run first
pixs=[250;250]; 
c=10; a=2;
G=d_maker(pixs(1),c,a);
Ga = fft2(G);
m=10; %History Window
lratconst=(1./vars)/2;
nums=X_count;
%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\///////////////////////////////////%
sim=10000;
UCL=[7300];
xcen=220; ycen=220;
fm=[2,3,4];
fs=[0.05,0.06,0.07,0.08];
fault=combvec(fm,fs);
H=fspecial('gaussian',4,2);
ARL=zeros(1,size(fault,2));
I=zeros(pixs(1),pixs(2));
xx1=zeros(pixs(1),pixs(2));
eps=2;

tic
for j=1:size(fault,2)
    counter=zeros(1,sim-m);
    X_sum=zeros(1580,m);
                  for picc=1:m
%                         I=iso_Exp(Ga,pixs(1));
                        cd stain
                 xx1=iso_Exp(Ga,pixs(1));
                 I=stainmake(xx1,fault(1,j),fault(2,j),eps,xcen,ycen,H);
                 cd ../
                    for p=1:1580
                          X_sum(p,picc)=sum(I(ROIs(:,:,p)));
                        end
                  end
     
    rat_tmp=zeros(1,sim-m);
    X_sum_new=zeros(1580,m);
    muhat=zeros(1580,1);
    rats=zeros(1580,1);
    tt=zeros(1580,1);

for pic= m+1:sim    
    flag=0;
    while flag == 0
%         I=iso_Exp(Ga,pixs(1));
[j, pic]
                 cd stain
                 xx1=iso_Exp(Ga,pixs(1));
                 I=stainmake(xx1,fault(1,j),fault(2,j),eps,xcen,ycen,H);
                 cd ../
        [tt]=evalROI_mex(ROIs,I); 
        X_sum_new(:,1:m-1)=X_sum(:,2:m);
        X_sum_new(:,m)=tt(:,1);
        X_sum=X_sum_new;
        [rat_max]= maxim_mex(m,X_sum,mus,nums,lratconst); 
        rat_tmp=(max(rat_max));
        counter(pic-m)=counter(pic-m)+1;
                   if rat_tmp > UCL
                    flag=1;
                   end
    
    end

  end
     ARL(1,j)=sum(counter)/(sim-m);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
