function X= iso_Exp(Ga,n);
Z = randn(n) + sqrt(-1)*randn(n);
X = real(fft2((Ga.^.5).*Z/n));
end