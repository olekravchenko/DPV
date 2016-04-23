function y = box(x)

x_dim = length(x);

y = zeros(1, x_dim);

y(1: ceil(0.2 * x_dim)) = 0;
y(ceil(0.2 * x_dim)+1 : ceil(0.4 * x_dim)) = 1;
y(ceil(0.4 * x_dim)+1 : end) = 0;