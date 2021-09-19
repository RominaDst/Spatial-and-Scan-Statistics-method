function [Ga]= iso_Exp_new(n,c,a);
t1 = [0:1/n:1-1/n]; t2 = t1;
for i=1:n 
for j=1:n
G(i,j)=exp(-c*(sqrt(min(abs(t1(1)-t1(i)), ...
1-abs(t1(1)-t1(i)))^2 + min(abs(t2(1)-t2(j)), ...
1-abs(t2(1)-t2(j)))^2)))^a;
end;
end;
Ga = fft2(G); 