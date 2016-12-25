%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                   Numerical solution of ABH model 
%                   (1)        U_t + F(U)_x = f,
%                              U = [nd vd]';
%                   F(U) = [nd .* vd; 0.5 * vd .^ 2 + taud_d * log(nd)]';
%                           f = [0 (a/(b+abs(vi)^3)-1)E - alpha_0 v_d]'
%                           with initial condition 
%                           nd(x,0) = 1e-3, vd(x,0) = 0
%                       by Lax-Friedrichs scheme 
%                   (2)         w = [ne E]';
%                           w_x = [-neE/ti 1-nd-ne]
%                           with initial condition
%                           ne(0) = 0.999, E(0) = 4e-3
%                    coded by Oleg Kravchenko 2013.01.29
%                    UPD1: 23/04/2016: Comments updated
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Refs:
% [1] K. Avinash, A. Bhattacharjee, and S. Hu, "Nonlinear Theory of Void 
%       Formation in Colloidal Plasmas," Phys. Rev. Lett., Vol. 90, 
%       No. 7, 075001-1-4, 2003.
%       DOI: 10.1103/PhysRevLett.90.075001
% [2] C.S. Ng, A. Bhattacharjee, S. Hu, Z.W. Ma, and K. Avinash, 
%       "Generalizations of a nonlinear fluid model for void formation
%       in dusty plasmas," Plasma Phys. Control. Fusion, Vol. 49, 
%       1583-1597, 2007.
%       DOI: http://dx.doi.org/10.1088/0741-3335/49/9/015
%       arXiv: http://arxiv.org/abs/1109.1039
% [3] O.V. Kravchenko, V.I. Pustovoit,"Numerical Simulation of Dynamics of 
%       Concentric Dusty Plasma Structures", Conference Presentation, 
%       WSWPA-2016, Moscow, Russia
%       https://www.researchgate.net/publication/301594302_Numerical_Simulation_of_Dynamics_of_Concentric_Dusty_Plasma_Structures
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all; clc;

% Set Path
addpath(genpath('D:\GitHub\DPV\ABH Model'));
% Set Default Figure Options
set(0,'DefaultTextInterpreter','latex')
set(0,'DefaultAxesFontSize',14)

%% Simulation Parameters
nx      = 2^6;                                  % Number of grid points
nL      = 6;                                    % Length of domain
dx      = nL/(nx-1);                            % Spatial step size
x       = 0:dx:nL;                              % Spatial grid
h2      = dx^2;                                 % Square of spatial step
cfl     = 0.45;                                 % CFL number for stability
sigma   = cfl;
dt      = sigma*dx;                             % Time step based on CFL number
tfin    = 90;                                   % Final time
t       = 0.0;                                  % Current time value
step    = 0;                                    % Current time step
oft     = 0.15;                                 % Figure offset parameter
%% Parameters in Model
td      = 1e-3;                                 % T_dT_i / T_e^2Z_d ratio
ti      = 0.125;                                % T_i/T_e ratio
cd      = sqrt(td);                             % Dust sound velocity
a       = 7.5;                                  % Constant in ion-drag force
b       = 1.6;                                  % Constant in ion-drag force
alpha0  = 2.0;                                  % Normalized dust-neutral collision frequency
mu      = 1.5;                                  % Constant mobility
%% Intial conditions
nd      = 1e-3*ones(1,nx);                      % Dust density
vd      = zeros(1,nx);                          % Dust velocity
vi      = zeros(1,nx);                          % Ion velocity
E       = zeros(1,nx);                          % Electrical field currency
U_star  = zeros(2,nx);                          % Vector function contains nd, vd 
f       = zeros(2,nx);                          % Right hand side function
U       = [nd; vd];
L = 1:nx-1; R = 2:nx;                           % Indicies arrays

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main Loop
% Set Figure
figure('color','w')
while t<tfin
    % Compute Delta U's
    % dU = U(R)-U(L) at the cell boundaries { x_{i},x_{i+1}, ... }
    dU = U(:,R)-U(:,L);
    % Define F and f vector (right side)
    F = [nd .* vd; 0.5 * vd .^ 2 + td * log(nd)];
    f(2,:) = (a ./ (b + abs(vi).^3) - 1).*E - alpha0*vd;    
    % Fluxes to the left and right of 'F_ {i+1/2}'
    FL = F(:,L); FR = F(:,R);
    
    nu = max([abs(vd+cd) abs(vd-cd)]);
    
    % Update time step
    t = t + dt;

    Flux = 0.5*(FL+FR - nu*dU);

    % Compute next time step
    for i = 2:nx-1
        U_star(:,i) = U(:,i) - sigma * (Flux(:,i) - Flux(:,i-1)) + dt*f(:,i);
    end   
    
    % BCs
    U_star(:,1)     = U(:,2);       % Neumann condition to the left
    U_star(:,nx)    = U(:,nx-1);     % Neumann condition to the right
    
    nd = U_star(1,:);
    vd = U_star(2,:);  
    
    % Compute ne, E functions
%     [ne, E] = runge4_solver(nL,nx,nd');
    [ne, E] = euler1_solver(nL,nx,nd',ti);
    E = E';
    vi = - mu * E;
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Recompute U vector according to w vector
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    f(2,:) = (a ./ (b + abs(vi).^3) - 1).*E - alpha0*vd;   
        
    % Update time step  
	dU = U_star(:,R)-U_star(:,L);
    Flux = 0.5*(FL+FR - nu*dU);
  
    % Compute next time step
    for i = 2:nx-1
        U(:,i) = U_star(:,i) - sigma * (Flux(:,i) - Flux(:,i-1)) + dt*f(:,i);
    end   
     
    % BCs
    U(:,1)	 = U(:,2);         % Neumann condition to the left
    U(:,end) = U(:,end-1);     % Neumann condition to the right
    
    nd = U(1,:);
    vd = U(2,:);
    
    % Limiting nd > 1 ---> nd = 1
    nd(nd>1) = 1;
      
    % New time step
    step = step + 1;    
    disp(['dt = ' num2str(t,'%2.2f\n') ' over ' num2str(tfin) '  step:' num2str(step)])
    
%     %% Plots
%     subplot(2,2,[3,4]); plot(x,nd,'k','LineWidth',1); 
%     xlabel('$x$'); ylabel('$n_d$'); title('$n_d(x,t$)');
%     xlim([0 nL]); ylim([-oft 1+oft])    
%     % Timer 1: current time step value
%     text('Interpreter','latex',...
%     'String',['time step: ' num2str(step) '/' num2str(length(dt:dt:tfin))],...
%     'Position',[0.2 .6],...
%     'FontSize',10, 'color','b', 'interpreter','latex')
%     % Timer 2:  current time value
%     text('Interpreter','latex',...
%     'String',['time: ' num2str(t,'%2.2f\n') ' : ' num2str(tfin)],...
%     'Position',[.6 .5],...
%     'FontSize',10, 'color','b', 'interpreter','latex')
%     
%     subplot(2,2,2); plot(x,vd,'k','LineWidth',1); 
%     xlabel('$x$'); ylabel('$v_d$'); title('$v_d(x,t)$');  
%     xlim([0 nL]); 
% 
% 	subplot(2,2,1); plot(x,ne,x,E); %,'k','LineWidth',1); 
%     xlabel('$x$'); ylabel('$n_e,\,E$'); title('$n_e(x,t),\,E(x,t)$');
%     xlim([0 nL]); 
%     
%     % Draw data on current Figure
%     drawnow      
end