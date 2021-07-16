clear;clc;
load ROI_Lin;
load mu0;
n = 256;
c=10; a=2;
G=d_maker(n,c,a);
Ga = fft2(G);
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
fc=[2,5,8,12];
fa=[0.25,0.5,1,1.5];
fault=combvec(fc,fa);

ARL=zeros(1,size(fault,2));
tmp=zeros(n,n);
xx=zeros(n,n);

for jj=1:size(fault,2)
    counter=zeros(1,nos);
     Ga=iso_Exp_new(n,fault(1,jj),fault(2,jj));
   for i=1:nos
       flag =0;
       xx1=zeros(n,n);
      
        while flag ==0
               Z = randn(n) + sqrt(-1)*randn(n);
               tmp = real(fft2((Ga.^.5).*Z/n));
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

