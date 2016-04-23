function [out1, out2] = runge4(a, b, N, alpha, nd)

%function rk4_systems() approximates the solutions of systems of m
%differential equations that are written in the form

%dy1/dt = f1(t,y1,y2,...,ym)
%dy2/dt = f2(t,y1,y2,...,ym)
%.
%.
%.
%dym/dt = fm(t,y1,y2,...,ym)

%with t in the interval [a; b] and the initial conditions are in the
%m-dimensional vector alpha
%as with function runge_kutta4(), the inputs are the endpoints a and b, the
%number of subdivisions N in the interval [a; b], and the initial
%conditions - but this time, the initial condition is a vector

%The algorithmic scheme in this file was drawn from the book of Burden & Faires
%Numerical Analysis, 7th Ed.

%Author: Oleg Kravchenko
%Date  : Mar. 2013

% result = 0;

m = size(alpha,1);
if m == 1
   alpha = alpha';
end

h = (b-a)/(N-1);        %the step size
t(1) = a;
w(:,1) = alpha;     %initial conditions

%%RGK4
% for i = 2:N
%    k1 = h*f(t(i-1), w(:,i-1), nd(i-1));
%    k2 = h*f(t(i-1)+h/2, w(:,i-1)+0.5*k1,nd(i-1));
%    k3 = h*f(t(i-1)+h/2, w(:,i-1)+0.5*k2,nd(i-1)); 
%    k4 = h*f(t(i-1)+h, w(:,i-1)+k3,nd(i-1));
%    w(:,i) = w(:,i-1) + (k1 + 2*k2 + 2*k3 + k4)/6;
%    t(i) = a + (i-1)*h;
% end

%%RGK2
for i = 2:N
   k1 = h*f(t(i-1), w(:,i-1), nd(i-1));
   k2 = h*f(t(i-1)+h, w(:,i-1)+0.5*k1, nd(i-1));
   w(:,i) = w(:,i-1) + 0.5 * (k1 + k2);
   t(i) = a + (i-1)*h;
end


out1 = t';
out2 = w';


%function relating the right-hand side of the differential equation
%it has to be changed accordingly to the problem at hand
%in this case, the system of differential equations is:
%dy1/dt = y2
%dy2/dt = -y1 - 2exp(t) + 1
%dy3/dt = -y1 - exp(t) + 1
%change it before proceeding to the command line
function dy = f(x, y, nd)
tmp = 1 - y(1) - nd;

if tmp < 0
    tmp = 0;
end
  dy = [ - 8 * y(1) * y(2); tmp];

   