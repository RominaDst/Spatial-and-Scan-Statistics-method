function [gam]=emprical_variogram(data0,resy0,resx0,p0,n0);

data0=double(data); %Image Data
data0_2=data0.^2; %Square of Image Data
ps0=numel(p0)/2; %Number of pairs
var_tmp0=zeros(ps0,1); %Average variogram for each column (pre-allocate)
for i=1:ps0
     var_tmp0(i)=sum(data0_2(:,p0(i,1))-2*data0(:,p0(i,1)).*data0(:,p0(i,2))+data0_2(:,p0(i,2)));
end

var_tmp0=var_tmp0/256;
var0=zeros(length(n0),1); %Variogram
n0=[0 n0]; %Number for each Distance (modified)
cn0=cumsum(n0); %Cummulative sum of number for each distance (used for indexing)
for i=2:length(n0)
     var0(i-1)=sum(var_tmp0((cn0(i-1)+1):cn0(i)))/n0(i);
end

gam=var0(1:resx0/2,1)/2;

end
% csvwrite('raw.vgm.csv',[np' hs' var0(1:resx0/2,1)/2 hor ver]);
