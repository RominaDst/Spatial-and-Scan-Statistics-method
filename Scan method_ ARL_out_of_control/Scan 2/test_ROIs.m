clc;clear;close all;
load('smallest_shift')
load('largest_shift')
pic2(pic2>0)=1;
load('ROIs');
for i=1:1580
%     i=1:numel(ROIs)/250/250
    ROIs2=double(ROIs(:,:,i));
    R=sum(sum(ROIs2));
    F=sum(sum(pic2));
    FR=sum(sum(ROIs2.*pic2));
    DSC(i)=2*FR/(F+R);
    imagesc(ROIs2);
    pause(0.1);
    clf;

end
    