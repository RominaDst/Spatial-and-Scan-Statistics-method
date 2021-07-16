function t=knot_maker(n,k)
%n is number of control points
%k is spline order (degree + 1)
    t=[];
    for i = 0:n+k
        if i >= 0 && i < k
            t=[t 0];
        elseif i >=  k && i <= n
            t=[t i-k+1];
        else
            t=[t n-k+2];
        end
    end
    t=t/max(t);
end
