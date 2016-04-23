function [ne, E] = newton_solver(n,h,tau_i,nd)

%% initial parameters
h2 = h^2;
ne = ones(1,n);    
ne(end) = 0;
E = zeros(1,n);
psiVector = zeros(1,n-2);
psiMatrixD = zeros(n-2,n-2);
epsilon = 1e-6;
r = 10*ones(1,n);

al = zeros(n-2,1);
am = zeros(n-2,1);
ar = zeros(n-2,1);

%% computation
iter_ind = 0;
while max(r) > epsilon   
    % psiVector definition
    psiVector(1:n-2) = ne(2:n-1) .* (ne(3:n) - 2*ne(2:n-1) + ne(1:n-2)) - ...
            0.25 * (ne(3:n) - ne(2:n-1)).^2 + ...
            (h2/tau_i) * ne(2:n-1).^2 .* (1 - nd(2:n-1) - ne(2:n-1));            

    % psiMatrixD Jacoby matrix definition
    for m = 1:n-2
        if m>1 && m<n-2
            psiMatrixD(m,m-1) = ne(m) + 0.5 * (ne(m+1) - ne(m-1));
            psiMatrixD(m,m+1) = ne(m) - 0.5 * (ne(m+1) - ne(m-1));
            psiMatrixD( m,m ) = ne(m+1) + ne(m-1) - 4*ne(m) + ... 
                2 * (h2/tau_i)*ne(m)*(1 - nd(m)) - 3 * (h2/tau_i)*(ne(m))^2;
        elseif m==1
            psiMatrixD(1,1) = ne(3) + ne(1) - 4*ne(2) + ... 
                2*ne(2)*(h2/tau_i) * (1 - nd(m)) - (3/tau_i)*h2*(ne(2))^2;
            psiMatrixD(1,2) = ne(2) - 0.5 * (ne(3) - ne(1));
        else 
            psiMatrixD(n-2,n-3) = ne(n-2) + 0.5 * (ne(n-1) - ne(n-3));
            psiMatrixD(n-2,n-2) = ne(n-1) + ne(n-3) - 4*ne(n-2) + ... 
                2*ne(n-2)*(h2/tau_i) * (1 - nd(m)) - (3/tau_i)*h2*(ne(n-2))^2;
        end; % if
    end % for m loop
    
    %% explicit solution
%     n_tmp(2:end-1) = ne(2:end-1) - (psiMatrixD^(-1) * psiVector')';

    %% tridiagonal solution
    am(1)       = ne(3) + ne(1) - 4*ne(2) + ... 
                2*ne(2)*(h2/tau_i) * (1 - nd(1)) - (3/tau_i)*h2*(ne(2))^2;
    am(2:n-3)   = ne(3:n-2) + ne(1:n-4) - 4*ne(2:n-3) + ... 
                2 * (h2/tau_i)*ne(2:n-3).*(1 - nd(2:n-3)) - 3 * (h2/tau_i)*(ne(2:n-3)).^2;
    am(n-2)     = ne(n-1) + ne(n-3) - 4*ne(n-2) + ... 
                2 * (h2/tau_i)*ne(n-2).*(1 - nd(n-2)) - 3 * (h2/tau_i)*(ne(n-2)).^2;
            
    al(1)       = 0;
    al(2:n-2)   = ne(2:n-2) + 0.5 * (ne(3:n-1) - ne(1:n-3));
    
    ar(1)       = ne(2) - 0.5 * (ne(3) - ne(1));
    ar(2:n-3)   = ne(2:n-3) - 0.5 * (ne(3:n-2) - ne(1:n-4));
    ar(n-2)     = 0;
    
    b = psiMatrixD * ne(2:n-1)' - psiVector';
    
    sol = tridiag(al,am,ar,b);  
    sol2 = [1 sol' 0]';  
    n_tmp = sol2;
    
    % error vector
    r = abs(n_tmp' - ne);   

    disp(['iteration: ' num2str(iter_ind) ' error Newton: ' num2str(max(r))])   

    ne(2:end-1) = n_tmp(2:end-1);
    iter_ind = iter_ind + 1;
end; % while


%% E field
for m=1:length(E)-1
    E(m+1) = E(m) + h * (1 - ne(m) - nd(m));
end

