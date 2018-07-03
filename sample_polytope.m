function X = sample_polytope(nsamp, Aeq, beq, lb, ub)
%SAMPLE_POLYTOPE samples a convex polytope
%   X = SAMPLE_POLYTOPE(NSAMP,AEQ,BEQ,LB,UB) samples the polytope (P) defined
%   by AEQ, BEQ, LB and UB, NSAMP times.
%   X is the output matrices. each column is a sample taken from P
%   NSAMP is the number of samples taken from P
%   AEQ is the equality coefficient matrix
%   BEQ is the equality constraint values
%   LB and UB are the lower and upper bounds.
TOL = 1e-10;

n = size(Aeq, 2);
% Find basis for null space
[~, R, E] = qr(Aeq, 0);
diagr = abs(diag(R));
r = find(diagr >= TOL * diagr(1), 1, 'last'); %rank estimation
x2_inds = sort(E(1:r));
A2 = Aeq(:, x2_inds);
x1_inds = setdiff(1:n, x2_inds);
A1 = Aeq(:, x1_inds);

Acp = [ eye(n - r);
      -eye(n - r);
      A2 \ A1;
      -A2 \ A1];
  
bcp = [ ub(x1_inds);
      -lb(x1_inds);
      A2 \ beq - lb(x2_inds)
      ub(x2_inds) - A2 \ beq ];
  
bcp = bcp(~all(abs(Acp) < TOL, 2));
Acp = Acp(~all(abs(Acp) < TOL, 2), :);

options = optimoptions('linprog', 'Display', 'off');
[~, ~, exitflag] = linprog(zeros(size(Acp, 2), 1), Acp, bcp, [], [], [], [], [], options);
if exitflag == -2
    fprintf('Warning: Initial bounds may be infeasible.\n');
end

cp_options.isotropic = 0;   
cp_options.method = 'gibbs';
x1mat = cprnd(nsamp, sparse(Acp), bcp, cp_options)';

X = zeros(n, nsamp);
X(x1_inds, :) = x1mat;
X(x2_inds, :) = A2 \ (repmat(beq, 1, nsamp) - A1 * x1mat);
end