function y = conf_int_fun(x, opt, emumodel, ldata, cdata, iparam, isul)

CONF_LEVEL = 0.95;
MAX_VAL = 1e7;

global mids;
global dmids;

nrxn = size(emumodel.S, 2);
nmetab = length(emumodel.metabs);
nsrc_iso = length(ldata.src_isos);
nsrc_iso_fixed = nnz(~isnan(ldata.src_fracs));
nlmetab = length(ldata.metabs);

mids = cell(length(emumodel.emus), 1);
dmids.v = cell(length(emumodel.emus), nrxn);
dmids.conc0 = cell(length(emumodel.emus), nmetab);
dmids.src_fracs = cell(length(emumodel.emus), nsrc_iso);

Aeq = [ emumodel.Afixed sparse(size(emumodel.Afixed, 1), nmetab + nsrc_iso + nlmetab);
    emumodel.S(cdata.inds, :) sparse(length(cdata.inds), nmetab + nsrc_iso + nlmetab);                             
    emumodel.S(setdiff(1:nmetab, cdata.inds), :) sparse(nmetab - length(cdata.inds), nmetab + nsrc_iso + nlmetab) 
    sparse(nsrc_iso_fixed, nrxn + nmetab) sparse(1:nsrc_iso_fixed, find(~isnan(ldata.src_fracs)), ones(nsrc_iso_fixed, 1), nsrc_iso_fixed, nsrc_iso) sparse(nsrc_iso_fixed, nlmetab)];
beq = [ emumodel.bfixed;
    zeros(length(cdata.inds), 1);
    zeros(nmetab - length(cdata.inds), 1);
    ldata.src_fracs(~isnan(ldata.src_fracs)) ];

lb = [ zeros(nrxn, 1);
    zeros(nmetab, 1);
    zeros(nsrc_iso, 1);
    zeros(nlmetab, 1) ];
ub = [ MAX_VAL * ones(nrxn, 1);
    MAX_VAL * ones(nmetab, 1);
    ones(nsrc_iso, 1);
    ones(nlmetab, 1) ];
lb(nrxn + cdata.inds) = cdata.x;
ub(nrxn + cdata.inds) = cdata.x;

xstar = zeros(nrxn + nmetab + nsrc_iso + nlmetab, 1);
xstar(1:length(opt.vnet)) = opt.vxch + max(opt.vnet, 0);
xstar(length(opt.vnet) + (1:nnz(emumodel.reverse))) ...
    = opt.vxch(emumodel.reverse) - min(opt.vnet(emumodel.reverse), 0);
xstar(nrxn + (1:nmetab)) = opt.conc0;
xstar(nrxn + nmetab + (1:nsrc_iso)) = opt.src_fracs;
xstar(nrxn + nmetab + nsrc_iso + (1:nlmetab)) = opt.dilutions;

fval0 = INST_sqerror(xstar, emumodel, ldata);
n = sum(cellfun(@(x) numel(x(2:end, :)), ldata.mids));
p = nrxn - rank(full(Aeq(:, 1:nrxn))) + nmetab + nsrc_iso + nlmetab;
thresh = fval0 + fval0 / (n - p) * finv(CONF_LEVEL, 1, n - p);
fprintf('thresh = %g\n', thresh);

if ~isul
    x = -x;
end

if x <= 0
    y = fval0 - thresh;
    return
end

if iparam <= length(opt.vnet)
    if emumodel.reverse(iparam)
        rev_inds = find(emumodel.reverse);
        Aeqp = [ Aeq;
            sparse([1 1], [iparam length(opt.vnet) + find(rev_inds == iparam)], [1 -1], 1, length(xstar)) ];
    else
        Aeqp = [ Aeq;
            sparse(1, iparam, 1, 1, length(xstar)) ];
    end
    beqp = [ beq;
        opt.vnet(iparam) + x ];
    A = [];
    b = [];
elseif iparam <= length(opt.vnet) + nnz(emumodel.reverse)
    % exchange flux; assume sign of net flux stays the same
    rev_inds = find(emumodel.reverse);
    if opt.vnet(rev_inds(iparam - length(opt.vnet))) > 0
        Aeqp = [ Aeq;
            sparse(1, iparam, 1, 1, length(xstar)) ];
        A = [ -sparse(1, rev_inds(iparam - length(opt.vnet)), 1, 1, length(xstar)) ];
    else
        Aeqp = [ Aeq;
            sparse(1, rev_inds(iparam - length(opt.vnet)), 1, 1, length(xstar)) ];
        A = [ -sparse(1, iparam, 1, 1, length(xstar)) ];
    end
    beqp = [ beq;
            opt.vxch(iparam - length(opt.vnet)) + x ];
    b = -(opt.vxch(iparam - length(opt.vnet)) + x);
end

qpopts = optimset('Display', 'off', 'Algorithm', 'interior-point-convex', 'TolCon', 1e-9);
[x0, ~, exitflag] = quadprog(eye(length(xstar)), -xstar, [], [], Aeqp, beqp, lb, ub, [], qpopts);
if exitflag == 1
    fval1 = INST_sqerror(x0, emumodel, ldata);
else
    fval1 = Inf;
end

if ~isinf(fval1)
    funcs.objective = @(x) INST_sqerror(x, emumodel, ldata);
    funcs.gradient = @(x) INST_sqerror_grad(x, emumodel, ldata);
    funcs.constraints = @(x) [Aeqp; A] * x;
    funcs.jacobian = @(x) [Aeqp; A];
    funcs.jacobianstructure = @(x) [Aeqp; A];
    options.lb = lb;
    options.ub = ub;
    options.cl = [ beqp; -Inf(size(b)) ];
    options.cu = [ beqp; b ];
    options.ipopt.hessian_approximation = 'limited-memory';
    [x, info] = ipopt(x0, funcs, options);
    fval = INST_sqerror(x, emumodel, ldata);
else
    fval = fval1;
end

y = fval - thresh;

end