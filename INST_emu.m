function INST_emu(analysis_name, nloop, init_fracs_zero, seed)
%INST_EMU find optimal values for the metabolic fluxes, metabolite
%concentrations, substrate labeling states as well as dilution factors.
%   INST_EMU(ANAYSIS_NAME, NLOOP, INIT_FRACS_ZERO, SEED) generates 
%   result files based on the network and the metabolic labeling data,
%   starting from a random initial point. The result files are saved in the
%   folder results\ANALYSIS_NAME
%   
%   ANALYSIS_NAME is the name of the folder including all required data
%   files.
%   NLOOP is the number of iterations
%   INIT_FRACS_ZERO is a binary value indicating the labeling state of
%   metabolites on time zero
%   SEED is the random generator seed, based on which the inital points
%   will be generated
%
%   A test folder is included when downloading INST-FLUX which can be used
%   as a template for new data sets.  testgit




NETWORK_FNAME = 'network.csv';
ATOM_SEP = '_';

TOL = 1e-10;

MIN_VAL = 1e-6;
MAX_VAL = 1e6;
MAX_V0_MULT = 2;

if nargin < 2
    nloop = 1;
end

if nargin < 3
    init_fracs_zero = true;
end

if nargin < 4;
    rng('shuffle');
else
    rng(seed);
end

disp(analysis_name);

ldata = read_ldata(analysis_name);
if isempty(ldata)
    return
end

bnd_metabs = unique(cellfun(@(x) x(1:(strfind(x, ATOM_SEP) - 1)), ldata.src_isos, 'UniformOutput', false));
emumodel = generate_emumodel([analysis_name '/' NETWORK_FNAME], ldata.emus, bnd_metabs);  

if isempty(emumodel.Afixed)
    fprintf('Error: There are no fixed fluxes in the network.\n');
    return
end

% Determine set of "non-zero" isotopomers, i.e. isotopomers whose initial
% fractional composition of the pool may be non-zero.
nz_isos = ldata.src_isos;
if ~init_fracs_zero
    n = regexp(emumodel.emus, '(?<metab>\S+)_(?<atom_nums>\S+)', 'names');
    n = n(cellfun('length', n) == 1);
        
    emu_metabs = cellfun(@(x) x.metab, n, 'UniformOutput', false);    
    emu_atom_nums = cellfun(@(x) x.atom_nums - '0', n, 'UniformOutput', false);
        
    for metab = emumodel.metabs
        atom_nums = emu_atom_nums(strcmp(metab, emu_metabs));
        
        atoms = false(1, max(cellfun(@(x) max(x), atom_nums)));
        for i = 1:length(atom_nums)
            atoms(atom_nums{i}) = true;
        end
        
        for x = 1:(2^nnz(atoms) - 1)
            iso_sub = zeros(1, length(atoms));
            iso_sub(atoms) = dec2binvec(x, nnz(atoms));
            iso = strcat(metab, ATOM_SEP, sprintf('%d', iso_sub(end:-1:1)));
            
            if ~ismember(iso, nz_isos)
                nz_isos(end + 1) = iso;
            end
        end
    end
end

% Read pool concentrations
cdata = read_conc(analysis_name, emumodel, bnd_metabs);  

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
conc_mat2 = -sparse(find(li), cellfun(@(x) find(strcmp(x, ldata.metabs)), cdata.metabs(li)), cdata.x(li), length(cdata.inds), nlmetab);

Aeq = [ emumodel.Afixed sparse(size(emumodel.Afixed, 1), nmetab + nnz_iso + nlmetab);
    emumodel.S sparse(nmetab, nmetab + nnz_iso + nlmetab);                             
    sparse(length(cdata.inds), nrxn) conc_mat1 sparse(length(cdata.inds), nnz_iso) conc_mat2;
    sparse(nnz_iso_fixed, nrxn + nmetab) sparse(1:nnz_iso_fixed, find(~isnan(ldata.src_fracs)), ones(nnz_iso_fixed, 1), nnz_iso_fixed, nnz_iso) sparse(nnz_iso_fixed, nlmetab)];
beq = [ emumodel.bfixed;
    zeros(nmetab, 1);
    cdata.x .* ~li;
    ldata.src_fracs(~isnan(ldata.src_fracs)) ];

lb = [ -TOL * ones(nrxn, 1);
    MIN_VAL * ones(nmetab, 1);
    zeros(nnz_iso, 1);
    zeros(nlmetab, 1) ];
ub = [ MAX_VAL * ones(nrxn, 1);
    MAX_VAL * ones(nmetab, 1);
    ones(nnz_iso, 1);
    ones(nlmetab, 1) ];

options = optimoptions('linprog', 'Display', 'off');
[v, ~, exitflag] = linprog(ones(nrxn, 1), [], [], [emumodel.Afixed; emumodel.S], [emumodel.bfixed; zeros(nmetab, 1)], zeros(nrxn, 1), [], [], options);
if exitflag == -2
    fprintf('Error: Stoichiometric constraints are infeasible.  Check the network.\n');
    return;
end

n = regexp(nz_isos, '(?<metab>\S+)_(\S+)', 'names');
nz_metabs = unique(cellfun(@(x) x.metab, n, 'UniformOutput', false));

total_frac = ones(length(nz_metabs), 1);
n = regexp(ldata.src_isos, '(?<metab>\S+)_(\S+)', 'names');
for i = 1:length(ldata.src_fracs)
    if ~isnan(ldata.src_fracs(i))
        tf = strcmp(n{i}.metab, nz_metabs);
        total_frac(tf) = total_frac(tf) - ldata.src_fracs(i);
    end
end

num_frac = zeros(length(nz_metabs), 1);
n = regexp(nz_isos, '(?<metab>\S+)_(\S+)', 'names');
for i = 1:length(n)
    if ~ismember(nz_isos(i), ldata.src_isos(~isnan(ldata.src_fracs)))
        tf = strcmp(n{i}.metab, nz_metabs);
        num_frac(tf) = num_frac(tf) + 1;
    end
end

A = sparse(length(nz_metabs), nrxn + nmetab + nnz_iso + nlmetab);
for imetab = 1:length(nz_metabs)
    A(imetab, nrxn + nmetab + (1:nnz_iso)) ...
        = strncmp(strcat(nz_metabs(imetab), ATOM_SEP), nz_isos, length(nz_metabs(imetab)) + 1)';
end
b = ones(length(nz_metabs), 1);

max_v0 = MAX_V0_MULT * max(v);
max_conc0 = max(cdata.x);

v0mat = sample_polytope(nloop, [ emumodel.Afixed; emumodel.S ], ...
    [ emumodel.bfixed; zeros(nmetab, 1) ], zeros(nrxn, 1), max_v0 * ones(nrxn, 1));
conc0mat = sample_polytope(nloop, conc_mat1, cdata.x, zeros(nmetab, 1), ...
    max_conc0 * ones(nmetab, 1));

x0 = zeros(nrxn + nmetab + nnz_iso + nlmetab, 1);
for iloop = 1:nloop
    % Generate random starting point in bounds and project onto feasible
    % polytope.   
    
    v0 = v0mat(:, iloop);
    conc0 = conc0mat(:, iloop);     
    
    fracs0 = zeros(length(nz_isos), 1);
    fracs0(find(~isnan(ldata.src_fracs))) = ldata.src_fracs(~isnan(ldata.src_fracs));
    for i = 1:length(total_frac)
        metab = nz_metabs{i};
        
        fracs = rand(num_frac(i), 1) * total_frac(i);
        fracs = sort(fracs);
        fracs = diff([0; fracs]);
        fracs0(strncmp(strcat(metab, ATOM_SEP), nz_isos, length(metab) + 1) ...
            & ~cellfun(@(x) ismember(x, ldata.src_isos(~isnan(ldata.src_fracs))), nz_isos)) ...
            = fracs;
    end
    
    x0(1:nrxn) = v0;
    x0(nrxn + (1:nmetab)) = conc0;
    x0(nrxn + nmetab + (1:nnz_iso)) = fracs0;
    x0(nrxn + nmetab + nnz_iso + (1:nlmetab)) = ones(nlmetab, 1);
    
    if any(x0 < lb) || any(x0 > ub) || any(abs(Aeq * x0 - beq) > TOL)
        fprintf('x0 does not satisfy constraints.\n');
        continue;
    end
    f = INST_sqerror(x0, emumodel, ldata, nz_isos);
    fprintf('RMS standardized error of initial point: %g\n', ...
        sqrt(INST_sqerror(x0, emumodel, ldata, nz_isos) ...
        / sum(cellfun(@numel, ldata.mids))));
    fprintf('fval of initial point = %g \n', INST_sqerror(x0, emumodel, ldata, nz_isos))

    tic;       
    options = optimset('Display', 'iter', 'GradObj', 'on', 'Algorithm', 'interior-point');
    [x, fval, exitflag, output] = fmincon(@(x) INST_sqerror(x, emumodel, ldata, nz_isos), ...
        x0, A, b, Aeq, beq, lb, ub, [], options);
    time_grad = toc;

    fprintf(['fval of final point = ' num2str(fval)])
    fprintf('RMS standardized error of final point: %g\n', ...
        sqrt(fval / sum(cellfun(@numel, ldata.mids))));

    if ~exist('results', 'dir')
        mkdir('results')
    end
    if ~exist(['results/' analysis_name], 'dir')
        mkdir('results', analysis_name);
    end
    save(['results/' analysis_name '/INST_emu_' datestr(now, 30)], ...
        'x', 'fval', 'exitflag', 'output', ...
        'emumodel', 'ldata', 'cdata', 'nz_isos');
end

end

