clc;clear;warning off all;fclose all;close all
%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\Parameters////////////////////////////////%
pixs=[250;250]; 
num_images_baseline=50; 
init_size=22; %Initial size (width and height) of serveillence region
grids=[10;10]; %Number of Section to scan...x and y directions
increment_size=4; %Increment size of serveillence region....increase in width and height
%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
ROI_maker(pixs,init_size,grids,increment_size);
load('ROIs')
%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ 
load('Sim_Baseline')
c=10; a=2;
G=d_maker(pixs(1),c,a);
Ga = fft2(G);
Sim_Baseline(pixs,num_images_baseline,ROIs,Ga);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Functions
function G=d_maker(n,c,a);
t1 = [0:1/n:1-1/n];
t2 = t1;
for i=1:n 
 for j=1:n
G(i,j)=exp(-c*(sqrt(min(abs(t1(1)-t1(i)), ...
1-abs(t1(1)-t1(i)))^2 + min(abs(t2(1)-t2(j)), ...
1-abs(t2(1)-t2(j)))^2)))^a;

 end;
end;

