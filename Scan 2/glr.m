function[rat_tmp]=glr(pic,m,sim,I,ROIs,nums,X_sum,mus,lratconst);
% for pic= m+1:sim
%   I=iso_Exp(Ga,pixs(1));
   muhat=zeros(1580,1);
    rats=zeros(1580,1);
    tt=zeros(1580,1);
%         for p=1:1580
%            tt(p,1)=sum(I(ROIs(:,:,p)));
%         end
%         
%          [tt]=evalROI_mex(ROIs,I);
         [tt]=evalROI(ROIs,I);
        
%   k= mod(pic,m);
%   if k==0
%       k=10;
%   end

   k= mod(pic-1,m)+1;
   X_sum(:,k)=tt(:,1);
       
    muhat=sum(X_sum,2)./(nums.*m );
    rats=lratconst.*nums.*(muhat-mus).^2;
%     rat_tmp(pic-m)=(max(rats));
    rat_tmp=(max(rats));
    
% end