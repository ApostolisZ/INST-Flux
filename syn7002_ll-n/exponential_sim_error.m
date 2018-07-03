<<<<<<< HEAD
function [exp_sq_error, dexp_sq_error] = exponential_sim_error(ab,t,mids,mids_std)
% syms a1 a2 a3 b1 b2 b3 t1 t2 t3 t4 t5 t6 t7 t8 t9 t10
% ab = [a1 a2 a3 b1 b2 b3];
% t = [t1 t2 t3 t4 t5 t6 t7 t8 t9 t10];
numvars = numel(ab);
a = ab(1:numvars/2);
b = ab(numvars/2+1:end);
nummeas = numel(t);
t = reshape(t,[],1);
b = reshape(b,1,[]);
M = t*b;
EM = exp(M);
a = reshape(a,[],1);
sim_mid = EM*a;
error = sim_mid' - mids;
error = reshape(error,[],1);
Sigma = inv(diag(mids_std));
exp_sq_error = error'*Sigma*error;

if nargout >1 
    dexp_sq_error(1:numvars/2,1) = sum(EM,1)';
    ta = t*a';
    db = sum(ta.*EM,1)';
    dexp_sq_error = [dexp_sq_error;db];
end
end

=======
function [exp_sq_error, dexp_sq_error] = exponential_sim_error(ab,t,mids,mids_std)
% syms a1 a2 a3 b1 b2 b3 t1 t2 t3 t4 t5 t6 t7 t8 t9 t10
% ab = [a1 a2 a3 b1 b2 b3];
% t = [t1 t2 t3 t4 t5 t6 t7 t8 t9 t10];
numvars = numel(ab);
a = ab(1:numvars/2);
b = ab(numvars/2+1:end);
nummeas = numel(t);
t = reshape(t,[],1);
b = reshape(b,1,[]);
M = t*b;
EM = exp(M);
a = reshape(a,[],1);
sim_mid = EM*a;
error = sim_mid' - mids;
error = reshape(error,[],1);
Sigma = inv(diag(mids_std));
exp_sq_error = error'*Sigma*error;

if nargout >1 
    dexp_sq_error(1:numvars/2,1) = sum(EM,1)';
    ta = t*a';
    db = sum(ta.*EM,1)';
    dexp_sq_error = [dexp_sq_error;db];
end
end

>>>>>>> 2024247f7ad30fed5e7c393c7eddc78b1a8690d7
