clc;clear;warning off all;fclose all;close all
%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\Parameters////////////////////////////////%
load('ROIs')
load('Sim_Baseline')
pixs=[250;250]; 
c=10; aa=2;
a=20; b=15;
G=d_maker(pixs(1),c,aa);
Ga = fft2(G);
[BlkCirc_row]=rho_maker(pixs(1),a,b);
m=10; %History Window
lratconst=(1./vars)/2;
nums=X_count;
%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\///////////////////////////////////%
eps=7.5;
xcen=128; ycen=128;
sim=10010;
UCL=[7300 ];
fm=[1,2,3];
fs=[0.1,0.15,0.2];
fault=combvec(fm,fs);

H=fspecial('gaussian',4,2);
ARL=zeros(1,size(fault,2));
I=zeros(pixs(1),pixs(2));
pic=zeros(pixs(1),pixs(2));
field1=zeros(pixs(1),pixs(2));

tic
for j=1:size(fault,2)
    counter=zeros(1,sim-m);
    X_sum=zeros(1580,m);
                  for picc=1:m
                 cd local_anis
                 field1= aniso_maker(pixs(1),BlkCirc_row);
                 pic=iso_Exp(Ga,pixs(1));
                 I=stainmake(pic,field1,fault(1,j),fault(2,j),eps,xcen,ycen,H);
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
    X_Sum_hist=X_sum;
for pics= m+1:sim    
    flag=0;
    [j,pics]
    while flag == 0
                 cd local_anis
                 field1= aniso_maker(pixs(1),BlkCirc_row);
                 pic=iso_Exp(Ga,pixs(1));
                 I=stainmake(pic,field1,fault(1,j),fault(2,j),eps,xcen,ycen,H);
                 cd ../
        [tt]=evalROI_mex(ROIs,I); 
        X_sum_new(:,1:m-1)=X_sum(:,2:m);
        X_sum_new(:,m)=tt(:,1);
        X_sum=X_sum_new;
        [rat_max]= maxim_mex(m,X_sum,mus,nums,lratconst); 
        rat_tmp=(max(rat_max));
        counter(pics-m)=counter(pics-m)+1;
                   if rat_tmp > UCL
                    flag=1;
                    X_sum=X_sum_hist;
                   end
    
    end

  end
     ARL(1,j)=sum(counter)/(sim-m)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
