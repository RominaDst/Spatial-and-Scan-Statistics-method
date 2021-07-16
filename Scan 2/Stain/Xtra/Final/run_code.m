function im=stain_maker(clc;clear;close all
pic=imread('pic.tif');
pic=double(rgb2gray(pic));
[y,x,z]=size(pic);
np=100;
p=3;
o=p+1;
nc=9;
spacer=1/(p+nc-2*p);
t=-p*spacer:spacer:(1+3*spacer);
eps=3;
pts=linspace(0,1,np);
flag=1;
while flag
    c=cpoints_maker(round(x/2),round(y/2),.05*min(x,y),nc+1,eps,o);
    stain=bspeval(p,c',t,pts);
    xres = linspace(1,x,x) ;
    yres = linspace(1,y,y) ;
    [X,Y] = meshgrid(xres,yres) ;

    [in,on] = inpolygon(X,Y,stain(1,:),stain(2,:)) ;
    im=in+on;
    im=~im;
    im=200*double(im);
    pic2=im+pic;

    figure()
    imagesc(pic2);colormap gray
    pause(1)
    close all
end


    
    
