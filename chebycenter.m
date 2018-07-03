function [c,r] = chebycenter(A,b)
%CHEBYCENTER Compute Chebyshev center of polytope Ax <= b.
%  The Chebyshev center of a polytope is the center of the largest
%  hypersphere enclosed by the polytope. 
%  Requires optimization toolbox.

[n,p] = size(A);
an = sqrt(sum(A.^2,2));
A1 = [ A an ];
f = zeros(p+1,1);
f(p+1) = -1;

options = optimset;
options = optimset(options,'Display', 'off');
c = linprog(f,A1,b,[],[],[],[],[],options);
r = c(p+1);
c = c(1:p);
end
