clc;clear;warning off all;fclose all;close all
%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\Parameters////////////////////////////////%
load('ROIs')
load('Sim_Baseline')
pixs=[250;250]; 
c=10; a=2;
m=10; %History Window
lratconst=(1./vars)/2;
nums=X_count;
%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\///////////////////////////////////%
sim=1010;
UCL=[7300 ];
fc=[2,5];
fa=[0.25,0.5,1];
fault=combvec(fc,fa);
ARL=zeros(1,size(fault,2));
Z=zeros(pixs(1),pixs(2));

for j=1:size(fault,2)
    counter=zeros(1,sim-m);
     Ga=iso_Exp_new(pixs(1),fault(1,j),fault(2,j));
    X_sum=zeros(1580,m);
                  for picc=1:m
%                         I=iso_Exp(Ga,pixs(1));
                Z = randn(pixs(1)) + sqrt(-1)*randn(pixs(1));
               I = real(fft2((Ga.^.5).*Z/pixs(1)));
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
               Z = randn(pixs(1)) + sqrt(-1)*randn(pixs(1));
               I = real(fft2((Ga.^.5).*Z/pixs(1)));
               
        [tt]=evalROI(ROIs,I); 
        X_sum_new(:,1:m-1)=X_sum(:,2:m);
        X_sum_new(:,m)=tt(:,1);
        X_sum=X_sum_new;
        [rat_max]= maxim(m,X_sum,mus,nums,lratconst); 
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
