function cpoints=cpoints_maker(x,y,r,nc,eps,ka)
cpoints=zeros(nc,2);
rads=linspace(0,2*pi,nc);
for i=1:(nc-ka)
    cpoints(i,:)=[x+r*cos(rads(i)) y+r*sin(rads(i))]+normrnd(0,eps,1,2);
end
cpoints(nc-ka+1:nc,:)=cpoints(1:ka,:);
