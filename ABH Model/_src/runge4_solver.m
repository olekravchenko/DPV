function [ne, E, x] = runge4_solver(b,n,nd)

%% initialization
a = 0;

%% computational
eta(1) = 4e-4;
[x,y] = runge4(a, b, n, [.9999; eta(1)], nd);

% clear x
disp('runge4 solver processing...')

ne  = y(:,1);
E   = y(:,2); % * (1.2 / max(y(:,2)));