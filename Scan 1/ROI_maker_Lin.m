clc;clear;warning off all;close all;
pixs=[256;256]; %Number of Pixels to Reshape Image
init_size=31; %Initial size (width and height) of serveillence region
grids=[8;8]; %Number of Section to scan...x and y directions
 increment_size=0;   %Increment size of serveillence region....increase in width and height
spaces=pixs./grids; %Calculates Max Surviellance Size   
I=zeros(pixs(1),pixs(2));
center=zeros(grids(1),grids(2),2);
min_dis=zeros(grids(1),grids(2));
max_steps=min_dis;
ROIs=logical(repmat(I,1,1,64));
ROItmp0=logical(I);
ROItmp1=ROItmp0;
ROItmp2=ROItmp1;
ROItmp3=ROItmp1;
ROItmp4=ROItmp1;
stepper_counter=0;
stepper2_counter=0;
for i = 1:grids(1) %Run through x locations for center(i,j} of serveillance box
    for j = 1:grids(2) %Run through y locations for center(i,j} of serveillance box
        center(i,j,:)=[spaces(1)/2+(i-1)*spaces(1);spaces(2)/2+(j-1)*spaces(2)]; %Current center(i,j} of serveillance box
        min_dis(i,j)=min([squeeze(center(i,j,:))-init_size/2-1;pixs-squeeze(center(i,j,:))-init_size/2-1]); %Find Nearest Edge...Maximum amount the box can increase
%         tempsteps=[init_size:increment_size:init_size+min_dis(i,j)*2]; % Size Increments of the box
          tempsteps=[init_size];
        max_steps(i,j)=length(tempsteps);
        steps{i,j}=tempsteps;
%         stepper2_counter=0;
        for stepper=1:max_steps(i,j)
            stepper_counter=stepper_counter+1;
            
            if stepper==1
                stepper2_counter=stepper2_counter+1;
                I_2_test1_reg{i,j,stepper}=[center(i,j,1)-steps{i,j}(stepper)/2,center(i,j,1)+steps{i,j}(stepper)/2,center(i,j,2)-steps{i,j}(stepper)/2,center(i,j,2)+steps{i,j}(stepper)/2];
                ROItmp0(I_2_test1_reg{i,j,stepper}(1):I_2_test1_reg{i,j,stepper}(2),I_2_test1_reg{i,j,stepper}(3):I_2_test1_reg{i,j,stepper}(4))=1; %Take Pixels out of Picture in box
                ROIs(:,:,stepper_counter)=ROItmp0;
            else
                stepper2_counter=stepper2_counter+1;
                I_2_test1_reg{i,j,stepper}=[center(i,j,1)-steps{i,j}(stepper)/2,center(i,j,1)-steps{i,j}(stepper)/2+increment_size/2-1,center(i,j,2)-steps{i,j}(stepper)/2,center(i,j,2)+steps{i,j}(stepper)/2];
                I_2_test2_reg{i,j,stepper}=[center(i,j,1)+steps{i,j}(stepper)/2-increment_size/2+1,center(i,j,1)+steps{i,j}(stepper)/2,center(i,j,2)-steps{i,j}(stepper)/2,center(i,j,2)+steps{i,j}(stepper)/2];
                I_2_test3_reg{i,j,stepper}=[center(i,j,1)-steps{i,j}(stepper)/2+increment_size/2,center(i,j,1)+steps{i,j}(stepper)/2-increment_size/2,center(i,j,2)-steps{i,j}(stepper)/2,center(i,j,2)-steps{i,j}(stepper)/2+increment_size/2-1];
                I_2_test4_reg{i,j,stepper}=[center(i,j,1)-steps{i,j}(stepper)/2+increment_size/2,center(i,j,1)+steps{i,j}(stepper)/2-increment_size/2,center(i,j,2)+steps{i,j}(stepper)/2-1,center(i,j,2)+steps{i,j}(stepper)/2+increment_size/2-2];
                ROItmp1(I_2_test1_reg{i,j,stepper}(1):I_2_test1_reg{i,j,stepper}(2),I_2_test1_reg{i,j,stepper}(3):I_2_test1_reg{i,j,stepper}(4))=1; %Pixels on Top of box
                ROItmp2(I_2_test2_reg{i,j,stepper}(1):I_2_test2_reg{i,j,stepper}(2),I_2_test2_reg{i,j,stepper}(3):I_2_test2_reg{i,j,stepper}(4))=1; %Pixels on Bottom of box
                ROItmp3(I_2_test3_reg{i,j,stepper}(1):I_2_test3_reg{i,j,stepper}(2),I_2_test3_reg{i,j,stepper}(3):I_2_test3_reg{i,j,stepper}(4))=1;  %Pixels on Left of box
                ROItmp4(I_2_test4_reg{i,j,stepper}(1):I_2_test4_reg{i,j,stepper}(2),I_2_test4_reg{i,j,stepper}(3):I_2_test4_reg{i,j,stepper}(4))=1;  %Pixels on Right of box
                ROIs(:,:,stepper_counter)=ROIs(:,:,stepper2_counter-1)+ROItmp1+ROItmp2+ROItmp3+ROItmp4;
                ROItmp1=logical(I);
                ROItmp2=ROItmp1;
                ROItmp3=ROItmp1;
                ROItmp4=ROItmp1;
            end
        end
        ROItmp0=logical(I);
    end
end
clearvars -except ROIs;
save('ROI_Lin','ROIs');