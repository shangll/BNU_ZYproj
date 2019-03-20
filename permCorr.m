function permCorr()
x = [1,5,3,2,5,7,2];
y = [2,3,1,1,4,3,3];

r = corrcoef(x,y);
r_value = r(1,2);

if r_value>=0
    t = .95;
else
    t = .05;
end

threshold_value = bsp(x,y,t);

if abs(r_value) > abs(threshold_value)
    r_value
    threshold_value
    ['accept it']
end
if (r_value) <= abs(threshold_value)
     r_value
     threshold_value
     ['reject it']
end

end

function threshold_value = bsp(x,y,t)

for k = 1:1000
    n = numel(x);
    data_s = [x,y];
    ns = randperm(numel(data_s));
    nx = ns(1:n);
    ny = ns(n+1:numel(data_s));
    data_xs = data_s(nx);
    data_ys = data_s(ny);
    rs = corrcoef(data_xs, data_ys);
    rs_value(k) = rs(1,2);
    
end

tempV = sort(rs_value);
thresholdN = floor(t*length(tempV));
threshold_value = tempV(thresholdN);

end
