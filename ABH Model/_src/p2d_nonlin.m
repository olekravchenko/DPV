% Dirichlet boundary conditions.  Uses a uniform mesh with (n+2)x(n+2) total
% points (i.e, n interior grid points).
% Input:
%        pfunc : the RHS of poisson equation (i.e. the Laplacian of u).
%        bfunc : the boundary function representing the Dirichlet B.C.
%          a,b : the interval defining the square
%            n : n+2 is the number of points in either direction of the mesh.
% Ouput:
%            u : the numerical solution of Poisson equation at the mesh points.
%          x,y : the uniform mesh.
%
% modified by Oleg Kravchenko  12/09/2013
function [output1, output2] = p2d_nonlin(n, h, ne, nd)

% h = (b-a)/(n+1);   % Mesh spacing
% tvect = a:h:b;
% [x,y] = meshgrid(tvect);   % Uniform mesh, including boundary points.

% Compute u on the boundary from the Dirichlet boundary condition
% ub = zeros(n,n);
idx = 2:n+1;
idy = 2:n+1;
% West and East boundaries need special attention
% ub(:,1) = ua(id,x(idx,1), y(idy,1));             % Left Boundary
% ub(:,n) = ua(id,x(idx,n+2), y(idy,n+2));         % Right Boundary
% Now the North and South boundaries
% ub(1,1:n) = ub(1,1:n) + ua(id,x(1,idx),y(1,idy)); 
% ub(n,1:n) = ub(n,1:n) + ua(id,x(n+2,idx),y(n+2,idy));


% u = ua(id,x(idx,idy),y(idx,idy));
u = zeros(n,n);
% u = ub;

% Create the D2x and D2y matrices
% Full matrix version.  This could be made much faster by using Matlab's
% sparse matrix functions (see "spdiags" for more details).
z = [-2;1;zeros(n-2,1)];
D2x = 1/h^2*kron(toeplitz(z,z),eye(n));
D2y = 1/h^2*kron(eye(n),toeplitz(z,z));

dunorm  = inf;
tol     = 1e-3;
iter    = 0;
% BC
% uUp     = ua(id,x(1,1:n+2),y(1,1:n+2)) + 0*ones(1, n+2) + 0*tvect;
% uDown   = ua(id,x(n+2,1:n+2),y(n+2,1:n+2)) + 0*tvect;
% uLeft   = ua(id,x(2:n+1,1),y(2:n+1,1)) + 0*tvect(2:n+1)';
% uRight  = ua(id,x(2:n+1,n+2),y(2:n+1,n+2)) + 0*ones(n, 1) + 0*tvect(2:n+1)';
uUp     = [ne(1) ne ne(n)];
uDown   = zeros(1, n+2) + ne(n);
uLeft   = ne.';
uRight  = zeros(n, 1) + ne(n);
% nd      = [BC(2,1) BC(2,:) BC(2,n)],n+2,1);
Emat    = eye(n+2);
while dunorm>tol
%     u = [ua(id,x(1,1:n+2),y(1,1:n+2)); ...
%             [ua(id,x(2:n+1,1),y(2:n+1,1)) u ua(id,x(2:n+1,n+2),y(2:n+1,n+2))]; ...
%                 ua(id,x(n+2,1:n+2),y(n+2,1:n+2))];
    u =         [uUp; ...
            [uLeft u uRight]; ...
                uDown];
 
    phiVector = 4*del2(u,h) - 8 * (exp(u) + nd - Emat);
    u = reshape(u(idx,idy),n*n,1);
%     u = reshape(u,n*n,1);
    phiVector = reshape(phiVector(idx,idy),n*n,1);
%     phiVector = reshape(phiVector,n*n,1);
    
    phiMatrix = D2x+D2y - diag(exp(u));
    b = phiMatrix * u - phiVector;      
    un = phiMatrix \ b;
    
    dunorm = norm(un - u,inf);
    % Convert f to a vector using column reordering
%     b = reshape(b,n*n,1);
    disp(['Iteration:  ' num2str(iter) ' Norm:   ',num2str(dunorm)]);
    
     u = un;
     u = reshape(u,n,n);
     iter = iter + 1;
end
% Solve the system

% Convert u from a column vector to a matrix to make it easier to work with
% for plotting.
u = reshape(un,n,n);

% Append on to u the boundary values from the Dirichlet condition.

% u = [[ua(id,x(1,1:n+2),y(1,1:n+2))];...
%      [[ua(id,x(2:n+1,1),y(2:n+1,1))] u ...
%       [ua(id,x(2:n+1,n+2),y(2:n+1,n+2))]];...
%      [ua(id,x(n+2,1:n+2),y(n+2,1:n+2))]];
 
%  h = x(1,2) - x(1,1);

u =         [uUp; ...
        [uLeft u uRight]; ...
            uDown];
        
% outputs
output1 = exp(u);
output2 = - 0.125 * 4 * del2(u,h);