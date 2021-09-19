function [tt]=evalROI(ROIs,I);
tt=zeros(1580,1);
for p=1:1580
           tt(p,1)=sum(I(ROIs(:,:,p)));
end
 