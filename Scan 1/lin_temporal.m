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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% UCL=[ 3.32 3.33 3.35 ];
UCL=[ 3.342 3.345 3.348 ];
nos=1000;
for jj=1:length(UCL)
    counter=zeros(1,nos);
   for i=1:nos
       UCL(jj)
       i
       flag =0;
       xx1=zeros(n,n);
        while flag ==0
              xx=iso_Exp(Ga,n);
              xx1=(1-w).*xx1 + w.*(xx);
              S=zeros(1,64);
              k=0;
                 for j=1:64
                    ss=0;
                     asdf=xx1(ROIs(:,:,j));
                      asdf1=reshape(asdf,32,32);
                        cen=[16 16];
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
                   if k > UCL(jj)  
                    flag=1;
                   end
        end
   end
   ARL(jj)=sum(counter)/nos;
    MRL(jj)=median(counter);
end
save( 'ARL_Lin_temporal3.mat', 'ARL','UCL','nos','MRL'); 