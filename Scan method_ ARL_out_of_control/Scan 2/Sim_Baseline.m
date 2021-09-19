function []=Sim_Baseline(pixs,num_images_baseline,ROIs,Ga)
X_sum=zeros(1580,num_images_baseline);
z=zeros(250,250,num_images_baseline);
for pic=1:num_images_baseline 
    I=iso_Exp(Ga,pixs(1));
    z(:,:,pic)=I;
    for p=1:1580
        X_sum(p,pic)=sum(I(ROIs(:,:,p)));
         X2_sum(p,pic)=sum((I(ROIs(:,:,p)).^2));
        if pic==1
        X_count(p,pic)= numel(I(ROIs(:,:,p)));
        end
    end
end
mus=sum( X_sum,2)./(X_count.*num_images_baseline );
vars=sum( X2_sum,2)./(X_count.*num_images_baseline )-mus.^2;
save('Sim_Baseline','mus','vars','X_count','z');

