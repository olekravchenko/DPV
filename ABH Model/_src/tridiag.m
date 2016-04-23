function u = tridiag(AL,AM,AR,r)
% Tridiagonal solver from the text Numerical Recipes 
% AL is the subdiagonal
% AM is the main diagonal
% AR is the superdiagonal
% r is the right hand side
% u is the solution to Au=r

n = length(r);

bet = AM(1);
u = zeros(n-1,1);
gam = u;
u(1) = r(1)/bet;

for j=2:n
   gam(j) = AR(j-1)/bet;
   bet = AM(j)-AL(j)*gam(j);
   if (bet==0) 
       disp('trinumrec failed'), pause, 
   end
   u(j) = (r(j)-AL(j)*u(j-1))/bet;
end

for j=n-1:-1:1
   u(j) = u(j) - gam(j+1)*u(j+1);
end



