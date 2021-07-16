clc;clear;warning off all;fclose all;close all
%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\Parameters////////////////////////////////%
load('ROIs')
pixs=[250;250]; 
num_images_baseline=50; 
init_size=22; %Initial size (width and height) of serveillence region
grids=[10;10]; %Number of Section to scan...x and y directions
increment_size=4;   %Increment size of serveillence region....increase in width and height
load('Sim_Baseline')
c=10; a=2;
G=d_maker(pixs(1),c,a);
Ga = fft2(G);
Sim_Baseline(pixs,num_images_baseline,ROIs,Ga);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

