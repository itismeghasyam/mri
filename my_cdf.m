function y =  my_cdf(x)

x = x(:);

x = sort(x, 'descend');


for i=1:length(x)
    y(i) = sum(x(1:i));
end

if y(end) > 0
y = y/y(end);
end
end