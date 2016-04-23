function [ne, E_star] = shooting_solver(b,n,nd)

%% initialization
a = 0;
h = (b - a) / (n-1);

E = zeros(1,n);
E_star = E;

%% computational
B = 1e-6;
eta(1) = - 1;
eta(2) = - 10;
[x,y] = runge4(a, b, n, [1; eta(1)], nd);
Phi(1) = abs(y(end,1) - B);
[x,y] = runge4(a, b, n, [1; eta(2)], nd);
Phi(2) = abs(y(end,1) - B);

k = 2;
disp('shooting solver processing...')
while Phi(k) > 1e-6
    k = k + 1;
    eta(k) = eta(k-1) - (eta(k-1) - eta(k-2)) / (Phi(k-1) - Phi(k-2)) * Phi(k-1);  
    [x, y] = runge4(a, b, n, [1; eta(k)], nd);
    Phi(k) = abs(y(end,1) - B);   
    
    disp(['iteration: ' num2str(k) '          error: ' num2str(Phi(k))])
end

ne = y(:,1);
for m = 2:length(E)-1
    E(m+1) = E(m-1) + 2*h * (1 - ne(m) - nd(m));
    E_star(m+1) = E_star(m) + h * (1 - ne(m) - nd(m));
end