% clc;clear;close all
 function pic2=stainmake(pic,fm,fs,eps,xcen,ycen,H);
% pic=imread('pic.tif');
% pic=double(rgb2gray(pic));
[y,x,z]=size(pic);
%  fm=3;
% fs=0.07;
np=100;
p=3;
o=p+1;
nc=18;
spacer=1/(p+nc-2*p);
t=-p*spacer:spacer:(1+3*spacer);
% eps=2;
pts=linspace(0,1,np);
% flag=1;
% while flag
    c=cpoints_maker(xcen,ycen,fs*min(x,y),nc+1,eps,o);
    stain=bspeval(p,c',t,pts);
    xres = linspace(1,x,x) ;
    yres = linspace(1,y,y) ;
    [X,Y] = meshgrid(xres,yres) ;
    [in,on] = inpolygon(X,Y,stain(1,:),stain(2,:)) ;
%         on=conv2(double(on),H);
%       on=on(1:256,1:256);

    im=in+on;
%     im=~im;
    im=fm*double(im);
    im=conv2(im,H);
     im=im(1:250,1:250);
    pic2=im+pic;
%     figure()
%       imagesc(pic2);colormap bone;colorbar
end
%     pause(1)
%     close all
% end


    
    
