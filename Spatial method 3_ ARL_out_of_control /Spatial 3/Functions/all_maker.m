function x = all_maker(s)
if s==1
%%%%%%%%%%%%%%%%
n = 256; % size of simulated surface n*n
c=10; a=2; % parameters of covarince function
G=covaraince_maker(n,c,a); 
% generating covaraince matrix
Ga = fft2(G);
Z = randn(n) + sqrt(-1)*randn(n);
X = real(fft2((Ga.^.5).*Z/n));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%induced local mean shift
xcen=120; ycen=120; % location of shift on the x and y axis
fm=[4]; % magnitude of the mean shift
fs=[0.09]; % shift size ratio (comapre to the original suface size n*n)
cd stain % The folder "stain" contains the NURBS codes for generating the random blot shaped shift
H=fspecial('gaussian',4,2); % create two-dimensional Guassian filter and will be applied on the border of the shift area
x=stainmake(X,fm,fs,eps,xcen,ycen,H); % function that induces the genearted shift on the in-control surface
cd ..
end

if s==2
 %%%%%%%%%%%%%%%   
n = 256; % size of simulated surface n*n
c=10; a=2; % parameters of covarince function
ti=pi/4; % anisotropic direction
R=[ cos(ti), -sin(ti); sin(ti),cos(ti)]; %transformantion matrix
b1=1;
b2=4;
b=b1/b2; % The smoothness of the correlation is controlled by the magnitude of separation distance
Ga=Ga_maker(n,R,a,c,b1,b2); % generating anisotropic covariance matrix
Z = randn(n) + sqrt(-1)*randn(n);
x = real(fft2((Ga.^.5).*Z/n));
end

if s==3
m = 256; % size of simulated surface n*n
c=10; a=2; % parameters of covarince function
fm=2;% magnitude of conditional mean shift
lambda2=-1;lambda1=4; % weights of adjacency matrix
cd GMRF % contain codes for generating GMRF
[idx,D,condk,mu,C,Qcond]=GMRF_input(m,lambda1,lambda2);%generate correlation function by using the Cholesky decomposition of the precision matrix 
yk=fm.*ones(m,1); 
mucond=mu(idx)+(Qcond\(-D(idx,condk)*(yk-mu(condk)))); % the condistional mean surface
x=MRF(m,mucond,C, idx, D, condk,yk); % generate the simulated surface
cd ..    
end

if s==4
 n = 256; % size of simulated surface n*n
c=10; a=2; % parameters of covarince function
G=covaraince_maker(n,c,a); % generating covaraince matrix
Ga = fft2(G);
Z = randn(n) + sqrt(-1)*randn(n);
x = real(fft2((Ga.^.5).*Z/n));
end

if s==5
 n = 256; % size of simulated surface n*n
c0=0; c=10; aa=2; % parameters of in-control covarince function
a=20; b=15; %  parameters of Anisotropic  covarince function
G=covaraince_maker(n,c,aa); % generating covaraince matrix
Ga = fft2(G);
cd local_anis % folder contain codes for generating blot shaped shift with NURBS 
[BlkCirc_row]=rho_maker(n,a,b); %generate block circulant matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%induced local mean shift
xcen=120; ycen=120; % location of shift on the x and y axis
fm=[4]; % magnitude of the mean shift
fs=[0.2]; % shift size ratio (comapre to the original suface size n*n)
eps0=7.5; % parameter to control the shape of the stain
H=fspecial('gaussian',4,2);% Guassian 2D Filter
field1= aniso_maker(n,BlkCirc_row); % % generate field with covariance given by block circulant matrix
pic=iso_Exp(Ga,n); % generate in-control surface
x=stainmake(pic,field1,fm,fs,eps0,xcen,ycen,H); % combine loacl shift with anistropic shift  
  cd ..  
end


end

%%%%%%%%%%%%%%%%%%%%Functions%%%%%%%%%%%%%%%%%%%%%%%%%

function [Ga]= Ga_maker(n,R,a,c,b1,b2);
t1 = [0:1/n:1-1/n]; t2 = t1;
for i=1:n
    for j=1:n
        dx= min(abs(t1(1)-t1(i)), 1-abs(t1(1)-t1(i)));
        dy=min(abs(t2(1)-t2(j)), 1-abs(t2(1)-t2(j)));
        dd=R*[dx;dy];
        d(i,j)=sqrt((dd(1)/b1)^2 +(dd(2)/b2)^2);
        G(i,j)=exp(-c*d(i,j))^a;
    end;
end;
Ga = fft2(G);
end

function [y1]= MRF(m,mucond,C, idx, D, condk,yk);
 y1=zeros(m,m);
 z=zeros(length(idx),1);
 x=zeros(length(idx),1);
z = randn(m^2-length(condk),1); % z is N(0,I)
x = mucond+C'\z;
y=zeros(m^2,1);
y(idx)=x;
y(condk)=yk;
y1=reshape(y,m,m);
end

function [idx,D,condk,mu,C,Qcond]= GMRF_input(m,d1,d2);
% clc;clear all
% m = 256; % Grid is m x m
% d1 = 1; % Central stencil weight
% d2 = -0.25; % Neighbour weight

% Elements conditioned on and their values
cond=[ones(m,1),(1:1:m)'];
condk=zeros(size(cond,1),1);

for i=1:size(cond,1)
    condk(i)=ij2k(cond(i,1),cond(i,2),m);
end

% yk=1.*ones(m,1);
nels = m*(5*m-4); % Number of non-zero elements on grid

% Pre-allocate memory to form the sparse precision matrix

a = zeros(1,nels);
b = zeros(1,nels);
c = zeros(1,nels);

% Compute the links and weights for precision matrix
k=0;
for i=1:m
    for j=1:m
        A = findneigh(i,j,m);
        nnb = size(A,1);
        for h=1:nnb
            a(k+h)= ij2k(i,j,m);
            b(k+h)= ij2k(A(h,1),A(h,2),m);
            if h==1
                c(k+h) = d1;
            else
                c(k+h) = d2;
            end
        end
        k = k+nnb;
    end  
end
% Construct the precision matrix

D = sparse(a,b,c,m^2,m^2);
% Construct mean vector

mu=sparse([],[],[],m^2,1,0);
% Obtain remaining indices
idx=(1:1:m^2)';
idx(condk)=[];
% Conditional Precision matrix
Qcond = D(idx,idx);

% mucond=mu(idx)+(Qcond\(-D(idx,condk)*(yk-mu(condk)))); 
C = chol(Qcond,'lower'); 
end

function G=covaraince_maker(n,c,a);
t1 = [0:1/n:1-1/n];
t2 = t1;
for i=1:n 
    for j=1:n
        G(i,j)=exp(-c*(sqrt(min(abs(t1(1)-t1(i)), ...
        1-abs(t1(1)-t1(i)))^2 + min(abs(t2(1)-t2(j)), ...
        1-abs(t2(1)-t2(j)))^2)))^a;
    end;
  end;
end

function [BlkCirc_row]=rho_maker(n,a,b);
m=n; 
% size of covariance matrix is m^2*n^2
tx=[0:n-1]; ty=[0:m-1]; % create grid for field

rho=@(x,y)((1-x^2/a^2-x*y/(b*a)-y^2/b^2)...
*exp(-(x^2/a^2+y^2/b^2)));

Rows=zeros(m,n); Cols=Rows;
for i=1:n
for j=1:m
Rows(j,i)=rho(tx(i)-tx(1),ty(j)-ty(1)); % rows of blocks
Cols(j,i)=rho(tx(1)-tx(i),ty(j)-ty(1)); % columns
end
end
% create the first row of the block circulant matrix
% with circulant blocks and store it as a matrix suitable for fft2
BlkCirc_row=[Rows, Cols(:,end:-1:2);
Cols(end:-1:2,:), Rows(end:-1:2,end:-1:2)];
end

function [field1]= aniso_maker(n,BlkCirc_row);
m=n;
lam=real(fft2(BlkCirc_row))/(2*m-1)/(2*n-1);
lam(lam(:)<0)=0; lam=sqrt(lam);
% generate field with covariance given by block circulant matrix
F=fft2(lam.*complex(randn(2*m-1,2*n-1),randn(2*m-1,2*n-1)));
F=F(1:m,1:n); % extract sub-block with desired covariance
field1=real(F);
end


