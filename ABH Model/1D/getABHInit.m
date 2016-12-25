function [nd, vd, vi, E, U_star, f] = getABHInit(id, nx)

switch id
    case '01'
        nd      = 1e-3*ones(1,nx);                      % Dust density
        vd      = zeros(1,nx);                          % Dust velocity
        vi      = zeros(1,nx);                          % Ion velocity
        E       = zeros(1,nx);                          % Electrical field currency
        U_star  = zeros(2,nx);                          % Vector function contains nd, vd
        f       = zeros(2,nx);                          % Right hand side function
    otherwise
        nd      = 1e-1*ones(1,nx);                      % Dust density
        vd      = 1e-3*ones(1,nx);                      % Dust velocity
        vi      = zeros(1,nx);                          % Ion velocity
        E       = zeros(1,nx);                          % Electrical field currency
        U_star  = zeros(2,nx);                          % Vector function contains nd, vd
        f       = zeros(2,nx);                          % Right hand side function
end        