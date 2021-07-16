clc;clear;
b1=csvread('wang_beta_historical.csv');
t1=csvread('wang_theta_historical.csv');
t1=t1(:,2:3);
t=mean(t1);
bm=mean(b1);
%%%%%%%%%%%%%%
n = 2^8;c=10;aa=2;
G=d_maker(n,c,aa);
Ga = fft2(G);
Z= iso_Exp(Ga,n);
a=20; b=15;
G=d_maker(n,c,aa);
Ga = fft2(G);
[BlkCirc_row]=rho_maker(n,a,b);
[resy0,resx0]=size(Z);
[p0,n0]=variogram_maker1(resx0);
np=n0(2:resx0/2+1);
hs=[1:resx0/2];
%%%%%%%%%%%%%%%
tmp=zeros(2,2);
S2=zeros(2,2);
l=length(t1);
for k=1:l-1
    tmp=(1/(2*(l-1)))*(t1(k+1,:)-t1(k,:))'*(t1(k+1,:)-t1(k,:));
    S2=S2+tmp;
end
iS2=inv(S2);
%%%%%%%%%%%%%%
r=16;
index=[1:(n/r):n];
x=index;
x=repmat(x,1,r);
y=index;
y=repmat(y,r,1);
y=reshape(y,1,numel(y));
B=[ones(r^2,1) x' y'];
%%%%%%%%%%%%%%%%%%%%%
options = optimoptions(@fmincon,'GradObj','on');
nos=10000;
eps=7.5;
xcen=128; ycen=128;
UCL1=[17];
UCL2=[22]; 
fm=[1.5,2,3];
fs=[0.1,0.15,0.2];
fault=combvec(fm,fs);
H=fspecial('gaussian',4,2);
ARL=zeros(1,size(fault,2));
field1=zeros(n,n);
Z=zeros(n,n);
pic=zeros(n,n);
%%%%%%%%%%%%%%%%%%%%%

tic
 for j=1:size(fault,2)
    counter=zeros(1,nos);
   for i=1:nos
       flag =0;
      
         while flag == 0
              z1=zeros(r,r);
              theta=zeros(1,3);
              beta=zeros(1,3);
              cd local_anis
                   field1= aniso_maker(n,BlkCirc_row);
                 pic=iso_Exp(Ga,n);
                 Z=stainmake(pic,field1,fault(1,j),fault(2,j),eps,xcen,ycen,H);
                 cd ../
%               Z= iso_Exp(Ga,n);  %image generating
              %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
              for jj=1:r
                z1(jj,:)=Z(((jj-1)*(n/r))+1,index); %downsampling
              end
              z=reshape(z1,r^2,1);
              %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
              [gam]=emprical_variogram_mex(Z,resy0,resx0,p0,n0);  %estimating parameters
              fun=@(x)fun2(x,np,gam,hs);
              theta= fmincon(fun,[.5,1,12],[],[],[],[],[0 0 0],[],[],options);
              sf = fit([x', y'],z,'poly11');
              beta= coeffvalues(sf);
              %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
              covert=wang_exp(r,theta(1),theta(2),theta(3),n);
              icovert=inv(covert);
              fisher1=B'*icovert*B;
               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
               T1=(beta-bm)*fisher1*(beta-bm)';
               T2=(theta(2:3)-t)*iS2*(theta(2:3)-t)';
               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
               counter(i)=counter(i)+1;
                    if T1 > UCL1 || T2> UCL2
                    flag=1;
                   end
        end
 end

   ARL(1,j)=sum(counter)/nos;
   
end
toc




    