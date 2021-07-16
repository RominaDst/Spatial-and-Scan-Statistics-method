clc;clear;warning off all;fclose all;close all
%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\Parameters////////////////////////////////%
load('ROIs')
load('Sim_Baseline')
pixs=[250;250]; 
c=10; a=2;
G=d_maker(pixs(1),c,a);
% Ga = fft2(G);
ti=pi/4;
R=[ cos(ti), -sin(ti); sin(ti),cos(ti)];
m=10; %History Window
lratconst=(1./vars)/2;
nums=X_count;
%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\///////////////////////////////////%
sim=1000;
UCL=[7300];
b=[1,1,1,1,1,1;2,3,4,5,10,20];
ARL=zeros(1,size(b,2));
Z=zeros(pixs(1),pixs(2));


for j=1:size(b,2)
    counter=zeros(1,sim-m);
     Ga=Ga_maker(pixs(1),R,a,c,b(1,j),b(2,j));
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
    [j pic]
    while flag == 0
    
%         I=iso_Exp(Ga,pixs(1));
               Z = randn(pixs(1)) + sqrt(-1)*randn(pixs(1));
               I = real(fft2((Ga.^.5).*Z/pixs(1)));
               
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
