clear;clc;
load ROI_Lin;
load mu0;
n = 256;
c=10; aa=2;
a=20; b=15;
G=d_maker(n,c,aa);
Ga = fft2(G);
[BlkCirc_row]=rho_maker(n,a,b);
l=0.2; %spatial
w=0.1; %temporal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sl=0;
for dd=1:16
    sl=sl+(1-l)^dd;
end
xx1=zeros(n,n);
cen=[16 16];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
UCL=[ 3.345 ];
nos=10000;
eps=7.5;
xcen=220; ycen=220;
fm=[1,2,3];
fs=[0.1,0.15,0.2];
fault=combvec(fm,fs);

H=fspecial('gaussian',4,2);
ARL=zeros(1,size(fault,2));
tmp=zeros(n,n);
xx=zeros(n,n);
pic=zeros(n,n);
field1=zeros(n,n);

tic
for jj=1:size(fault,2)
    counter=zeros(1,nos);
   for i=1:nos
       flag =0;
       xx1=zeros(n,n);
        while flag ==0
         
            cd local_anis
%                  field1= aniso_exp(n,a,b);
                   field1= aniso_maker(n,BlkCirc_row);
                 pic=iso_Exp(Ga,n);
                 tmp=stainmake(pic,field1,fault(1,jj),fault(2,jj),eps,xcen,ycen,H);
                 cd ../
              xx1=(1-w).*xx1 + w.*(tmp);
              S=zeros(1,64);
              k=0;
                 for j=1:64
                    ss=0;
                     asdf=xx1(ROIs(:,:,j));
                      asdf1=reshape(asdf,32,32);
                              for d=1:16 
                                x=asdf1(cen(1)-d+1:cen(1)+d,cen(2)-d+1:cen(2)+d);
                                asdf1(cen(1)-d+1:cen(1)+d,cen(1)-d+1:cen(1)+d)=0;
                                mu1=mean(mean(nonzeros(x)));
                                ss=(((1-l)^d)/sl)*(mu1-(mu0(j,:)))+ss;
                              end
                              S(j)=ss;
                           end  
                         k=max(S)/sqrt(var(S));
               counter(i)=counter(i)+1;
                   if k > UCL  
                    flag=1;
                   end
        end
   end
   ARL(jj)=sum(counter)/nos;
end

