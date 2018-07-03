function calc_ci_samp(analysis_name, result_name, nloop)
%calc_ci_samp calculates  NLOOP minima based on data with added noise based
%on the standard deviation of the experimental measurements. The starting point for 
%the optimization problems is the optimal value, around which the
%confidence intervals need to be calculated. 

ATOM_SEP = '_';
TOL = 1e-10;
MIN_VAL = 1e-6;
MAX_VAL = 1e6;

if nargin < 3
    nloop = 1;
end

rng('shuffle');

results_dir = [ 'results/' analysis_name ];
load([results_dir '/' result_name]);

ldata0 = ldata; % bkup original ldata
cdata0 = cdata;

%
% Set up basic constants for equations
%

nrxn = size(emumodel.S, 2);
nmetab = length(emumodel.metabs);
nnz_iso = length(nz_isos);
nnz_iso_fixed = nnz(~isnan(ldata.src_fracs));
nlmetab = length(ldata.metabs);

temp_mat = cell(size(cdata.inds));
for i = 1:length(temp_mat)
    temp_mat{i} = ones(size(cdata.inds{i})) * i;
end
conc_mat1 = sparse(cell2mat(temp_mat), cell2mat(cdata.inds), 1, length(cdata.inds), nmetab);
li = cellfun(@(x) ismember(x, ldata.metabs), cdata.metabs);

lb = [ -TOL * ones(nrxn, 1);
    MIN_VAL * ones(nmetab, 1);
    zeros(nnz_iso, 1);
    zeros(nlmetab, 1) ];
ub = [ MAX_VAL * ones(nrxn, 1);
    MAX_VAL * ones(nmetab, 1);
    ones(nnz_iso, 1);
    ones(nlmetab, 1) ];

n = regexp(nz_isos, '(?<metab>\S+)_(\S+)', 'names');
nz_metabs = unique(cellfun(@(x) x.metab, n, 'UniformOutput', false));

A = sparse(length(nz_metabs), nrxn + nmetab + nnz_iso + nlmetab);
for imetab = 1:length(nz_metabs)
    A(imetab, nrxn + nmetab + (1:nnz_iso)) ...
        = strncmp(strcat(nz_metabs(imetab), ATOM_SEP), nz_isos, length(nz_metabs(imetab)) + 1)';
end
b = ones(length(nz_metabs), 1);

options = optimset('Display', 'iter', 'GradObj', 'on', 'Algorithm', 'interior-point');

for iloop = 1:nloop
    %
    %   Adding noise into ldata.mids based on ldata.mids_std
    %
    ldata = ldata0;
    for iion = 1:length(ldata.ions)
        ldata.mids{iion} = ldata.mids{iion} ...
            + randn(size(ldata.mids{iion})) .* ldata.mids_std{iion};
    end
    
    cdata = cdata0;
    cdata.x = cdata.x + randn(size(cdata.x)) .* cdata.x_std;
    
    conc_mat2 = -sparse(find(li), cellfun(@(x) find(strcmp(x, ldata.metabs)), cdata.metabs(li)), cdata.x(li), length(cdata.inds), nlmetab);
    
    Aeq = [ emumodel.Afixed sparse(size(emumodel.Afixed, 1), nmetab + nnz_iso + nlmetab);
        emumodel.S sparse(nmetab, nmetab + nnz_iso + nlmetab);
        sparse(length(cdata.inds), nrxn) conc_mat1 sparse(length(cdata.inds), nnz_iso) conc_mat2;
        sparse(nnz_iso_fixed, nrxn + nmetab) sparse(1:nnz_iso_fixed, find(~isnan(ldata.src_fracs)), ones(nnz_iso_fixed, 1), nnz_iso_fixed, nnz_iso) sparse(nnz_iso_fixed, nlmetab)];
    beq = [ emumodel.bfixed;
        zeros(nmetab, 1);
        cdata.x .* ~li;
        ldata.src_fracs(~isnan(ldata.src_fracs)) ];
    
    max_conc0 = max(cdata.x);
    conc0 = sample_polytope(1, conc_mat1, cdata.x, zeros(nmetab, 1), ...
        max_conc0 * ones(nmetab, 1));
    
    x0 = x;
    x0(nrxn + (1:nmetab)) = conc0;
    
    [x, fval, exitflag, output] = fmincon(@(x) INST_sqerror(x, emumodel, ldata, nz_isos), ...
        x0, A, b, Aeq, beq, lb, ub, [], options);
    
    if ~exist([results_dir '/' result_name], 'dir')
        mkdir([results_dir '/'], result_name);
    end
    save([results_dir '/' result_name '/ci_' datestr(now, 30)], ...
        'x', 'fval', 'exitflag', 'output');
end
end