clear;clc;
load ROI_Lin;
load mu0;
n = 256;
m=n;
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
UCL=[3.345];
nos=10000;
f=[1.5,2];
d=[4 2;-1 -0.5];
fault=combvec(f,d);
H=fspecial('gaussian',4,2);
ARL=zeros(1,size(fault,2));
tmp=zeros(n,n);
xx=zeros(n,n);

for jj=1:size(fault,2)
    cd GMRF
    [idx,D,condk,mu,C,Qcond]=GMRF_input(m,fault(2,jj),fault(3,jj));
    yk=fault(1,jj).*ones(m,1);
    mucond=mu(idx)+(Qcond\(-D(idx,condk)*(yk-mu(condk))));
    cd ../
    counter=zeros(1,nos);
   for i=1:nos
       [jj,i]
       flag =0;
       xx1=zeros(n,n);
        while flag ==0
         
            cd GMRF
                 tmp=MRF(m,mucond,C, idx, D, condk,yk);
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
