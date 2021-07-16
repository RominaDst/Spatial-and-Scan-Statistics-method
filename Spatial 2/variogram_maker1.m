function [p,n]=variogram_maker1(res)
resm1=res-1;
p=zeros(res*resm1/2,2);
counter=0;
n=resm1:-1:1;
for x1=1:resm1
    for x2=1:(res-x1)
        counter=counter+1;
        p(counter,:)=[x2 x2+x1];
    end
end



