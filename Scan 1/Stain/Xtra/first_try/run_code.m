clc;clear;close all
pic=imread('pic.tif');
[x,y]=size(pic);
n=100;
k=6;
nc=9; %Number of CPs -1....Closed Actual # + 1
t=knot_maker(nc,k);
% t=[0 0 linspace(0,1,10) 1 1];
eps=4;
ka=1;
c=cpoints_maker(round(x/2),round(y/2),.1*min(x,y),nc+1,eps,ka);
pts=linspace(0,1,n);
% for i=1:n
%     for j=1:(k-1)
%         left=find(t>pts(i),1);
%         left=left-1;
%         w=(pts(i)-t(left))/(t(left+j-1)-t(left));
%         Bleft=w*Bleft+(1-w)*Bright
%     end
% end
stain=bspeval(k-1,c',t,pts);
plot(c(:,1),c(:,2),'*');
hold on
plot(stain(1,:),stain(2,:))
axis equal
