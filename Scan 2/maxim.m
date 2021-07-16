function [ rat_max]= maxim(m,X_sum,mus,nums,lratconst);
rat_max=zeros(m,1);
for j=1:m
          muhat=sum(X_sum(:,m-j+1:m),2)./(nums.*j);
         rats=lratconst.*nums.*j.*(muhat-mus).^2;
          rat_max(j,1)=(max(rats));
end
end