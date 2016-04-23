function y = step(x)

x_dim = length(x);

y = zeros(1, x_dim);

y(1: ceil(.5 * x_dim)) = 0.01;
y(ceil(.5 * x_dim)+1 : end) = 1;