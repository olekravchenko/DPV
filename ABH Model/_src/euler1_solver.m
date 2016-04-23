function [ne, E, x] = euler1_solver(Lx,nx,nd,ti)

%% Initialization
x   = linspace(0,Lx,nx);    dx = Lx / (nx-1);
E   = 0*x';
ne  = E;
nd0 = 1e-3;
ne0 = 1 - nd0;
E0  = 4e-3;
ti_inv = 1/ti;

% IC
ne(1) = ne0;
E(1)  = E0;

% indexes
L = 1:nx-1;

for ind = L
    ne(ind+1)   =   ne(ind) - ti_inv * dx * ne(ind) .* E(ind);
    E(ind+1)    =   E(ind) + dx * (1 - ne(ind) - nd(ind));
end
