function [ne, E] = triblock_solver(L,n,tau_i,nd)

%% initial parameters
h = L / (n-1);
x = 0:h:L;

E = 1.2* x / L;
ne = 1 - x / L;
w = zeros(2*n,1);
w(1:2:end-1) = ne;
w(2:2:end) = E;
%%
epsilon = 1e-1;
r = 10*ones(1,n);
% G matrix, phi vector
G = zeros(2*(n-2),2*(n-2));
phi = zeros(2*(n-2),1);

%% computation
iter_ind = 0;
while max(r) > epsilon 
% while iter_ind <= 100 
    % G matrix
    for m = 2:2:2*(n-2)
        if m == 2
            G(1:2,1:2) = -2 * eye(2);
            G(1:2,3:4) = eye(2) - 0.5*h*[-E(2)/tau_i -ne(2)/tau_i; -1 0];
        elseif m == 2*(n-2)
            G(2*(n-2)-1:2*(n-2),2*(n-2)-3:2*(n-2)-2) = eye(2) + 0.5*h*[-E(n-1)/tau_i -ne(n-1)/tau_i; -1 0];        
            G(2*(n-2)-1:2*(n-2),2*(n-2)-1:2*(n-2)) = -2 * eye(2);        
        else
            ind = ceil(0.5*m);
            G(m-1:m,m-3:m-2) = eye(2) + 0.5*h*[-E(ind)/tau_i -ne(ind)/tau_i; -1 0];
            G(m-1:m,m-1:m) = -2 * eye(2);        
            G(m-1:m,m+1:m+2) = eye(2) - 0.5*h*[-E(ind)/tau_i -ne(ind)/tau_i; -1 0];
        end
    end % for G matrix
    % phi vector
    for m = 1:2:2*(n-2)
        phi(m) = 0;
        ind = ceil(0.5*m)+1;
        phi(m+1) = 0.5 * (-nd(ind-1) + nd(ind+1)) * h;
    end % for phi vector

    phi(1:2) = - (eye(2) + 0.5*h*[-E(2)/tau_i -ne(2)/tau_i; -1 0]) * [1; 0];
    phi(end-1:end) = - (eye(2) - 0.5*h*[-E(end-1)/tau_i -ne(end-1)/tau_i; -1 0]) * [0; 1.2];

%     phi = G * w(3:2*n-2) - phi;
    
    sol = triblocksolve(G,phi,n-2);

    %% output solution
    
    % error vector
    r = abs(w(3:2*n-2) - sol);   

    disp(['iteration: ' num2str(iter_ind) ' error Newton: ' num2str(max(r))])   
    
    w(3:2*n-2) = sol;
    ne(2:end-1) = sol(1:2:end-1);
    E(2:end-1) = sol(2:2:end);
    
    iter_ind = iter_ind + 1;
end; % while

%% output solution
ne(2:end-1) = sol(1:2:end-1);
E(2:end-1) = sol(2:2:end);