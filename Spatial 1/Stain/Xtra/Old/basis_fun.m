function N = basis_fun(d,n,k,u)
    nu = numel(u);
    N = zeros(nu,n+1);                               

    for col=1:nu                                    
        s = findspan(n, d, u(col), k);
        x=3;
        N(col,s+1-d:s+1) = basisfun_lee(s,u(col),d,k);
    end
end