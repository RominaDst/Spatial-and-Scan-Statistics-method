Rpath='C:\Program Files\R\R-3.4.1\bin';
RscriptFileName='C:\directory\to\fit-vario.R';% R code for fitting the variogram
sim=1;
n = 2^8;%image size
c=10;a=2;% in-control surface paramaters
theta=zeros(sim,4);
in_control_data=zeros(n^2,sim);
for i=1:sim
  Ga= iso_Exp_new(n,c,a); 
  X= iso_Exp(Ga,n);
  X=xx(:,i);
  X=reshape(X,n,n);
  emprical_variogram(X);% calulate the variogram
  RunRcode(RscriptFileName, Rpath);% fit the varigram and estimate the parameters theta(sill and range)
  theta(i,:)=csvread('theta.csv');
  in_control_data(:,i)=reshape(X,n^2,1);
end
 theta1=mean(theta);
 csvwrite('in_control_data.csv')