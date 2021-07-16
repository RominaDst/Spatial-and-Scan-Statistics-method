function nurbs = nrbmak_lee(coefs,knots) 
    nurbs.form   = 'B-NURBS'; 
    nurbs.dim    = 4; 
    np = size(coefs); 
    dim = np(1); 
    nurbs.number = np(2:3); 
    nurbs.coefs = repmat([0.0 0.0 0.0 1.0]',[1 np(2:3)]); 
    nurbs.coefs(1:dim,:,:) = coefs;   
    uorder = size(knots{1},2)-np(2); 
    vorder = size(knots{2},2)-np(3); 
    uknots = sort(knots{1}); 
    vknots = sort(knots{2}); 
    uknots = (uknots-uknots(1))/(uknots(end)-uknots(1)); 
    vknots = (vknots-vknots(1))/(vknots(end)-vknots(1)); 
    nurbs.knots = {uknots vknots}; 
    nurbs.order = [uorder vorder];
end
 
 
