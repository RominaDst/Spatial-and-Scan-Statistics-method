clc;clear;warning off all;fclose all;close all
%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\Parameters////////////////////////////////%
load('ROIs')
load('Sim_Baseline')
pixs=[250;250]; 
mm=pixs(1);
c=10; a=2;
G=d_maker(pixs(1),c,a);
Ga = fft2(G);
m=10; %History Window
lratconst=(1./vars)/2;
nums=X_count;
%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\///////////////////////////////////%
sim=10010;
UCL=[7300 ];
f=[1,1.5,2];
d=[4 2;-1 -0.5];
fault=combvec(f,d);
 
ARL=zeros(1,size(fault,2));
I=zeros(pixs(1),pixs(2));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
for j=1:size(fault,2)
    cd GMRF
    [idx,D,condk,mu,C,Qcond]=GMRF_input(mm,fault(2,j),fault(3,j));
    yk=fault(1,j).*ones(mm,1);
    mucond=mu(idx)+(Qcond\(-D(idx,condk)*(yk-mu(condk))));
    cd ../
    counter=zeros(1,sim-m);
    X_sum=zeros(1580,m);
                  for picc=1:m
%                         I=iso_Exp(Ga,pixs(1));
                       cd GMRF
                          I=MRF(mm,mucond,C, idx, D, condk,yk);
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
    [j, pic]
    while flag == 0
%         I=iso_Exp(Ga,pixs(1));
                 cd GMRF
                 I=MRF(mm,mucond,C, idx, D, condk,yk);
                 cd ../
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
toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
