function[rat_tmp]=glr(pic,m,sim,I,ROIs,nums,X_sum,mus,lratconst);
   muhat=zeros(1580,1);
    rats=zeros(1580,1);
    tt=zeros(1580,1);
    [tt]=evalROI(ROIs,I);
   k= mod(pic-1,m)+1;
   X_sum(:,k)=tt(:,1);      
    muhat=sum(X_sum,2)./(nums.*m );
    rats=lratconst.*nums.*(muhat-mus).^2;
    rat_tmp=(max(rats));
    