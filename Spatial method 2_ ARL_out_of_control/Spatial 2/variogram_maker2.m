function [data2,numnonzero]=variogram_maker2(data)
[res1,res2]=size(data);
res1p1=res1+1;
data2=NaN(res1+res2-1,min(res1,res2));
[res1a,res2a]=size(data2);
res1ap1=res1a+1;
numnonzero=[res1a:-2:1];

for i=1:res1
    data2(i:res1ap1-i,i)=[data(1:res1p1-i,i);data(res1p1-i,i+1:res2)'];
end


