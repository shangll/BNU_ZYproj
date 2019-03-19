function permCorr()
x = [66.5,58.5,69.5,77,65,68,71,60,72,71.5];
y = [.46,.43,.41,.45,.45,.45,.46,.47,.49,.44];
data_xy = [x(:),y(:)];

p1 = [2,1,2,1,2,1,2,2,1,1];
p2 = [9.31,9.43,8.84,11.62,10.20,10.45,10.31,9.50,10.80,10.82];
p = [p1(:),p2(:)];

[rho,pval] = partialcorr(data_xy,p);
r_value = rho(1,2);

threshold_value = bsp(x,y,p);

if abs(r_value) > abs(threshold_value)
    r_value
    threshold_value
    ['accept it']
end
if abs(r_value) <= abs(threshold_value)
    r_value
    threshold_value
     ['reject it']
end

end

function threshold_value = bsp(x,y,p)

for k = 1:1000
    n = numel(x);
    data_s = [x,y];
    ns = randperm(numel(data_s));
    nx = ns(1:n);
    ny = ns(n+1:numel(data_s));
    data_xs = data_s(nx);
    data_ys = data_s(ny);
    data_xy_s = [data_xs(:),data_ys(:)];
    
    [rho_s,pval_s] = partialcorr(data_xy_s,p);
    rs_value(k) = rho_s(1,2);
    
end

tempV = sort(rs_value);
thresholdN = floor(.95*length(tempV));
threshold_value = tempV(thresholdN);

end