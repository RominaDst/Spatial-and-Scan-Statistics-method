function summ=fitvaro(theta,np,gam,hs)
summ=0;
tmp=0;
    for i=1:length(np)
        var(i)= theta(1)+ theta(2)*(1-exp(-hs(i)/theta(3)));

        tmp=(np(i)/(var(i)^2))*((gam(i)-var(i))^2);
        summ=summ+tmp;
    end
end