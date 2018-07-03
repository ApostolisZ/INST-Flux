function [opt, emumodel, ldata, cdata, xbulk] = process_results(analysis_name, figures, write, estimate_conf)
%PROCESS_RESULTS processes calculated results produced by INST-EMU to find
%the best fit
%   [OPT, EMUMODEL,LDATA,CDATA,XBULK] = PROCESS_RESULTS(ANALYSIS_NAME,FIGURES,WRITE,ESTIMATE_CONF)
%   analyzes the minima calculated for ANALYSIS_NAME, selects and analyzes
%   the optimal result. 
%   ANALYSIS_NAME is the name of the folder containing the analysis data
%   FIGURES is a boolean operator determining whether figures will be
%   produced. Default value is TRUE
%   WRITE is a boolean operator determining whether the results will be
%   written in an output xlsx file in the folder ANALYSIS_NAME. Default
%   value is FALSE
%   ESTIMATE_CONF is a boolean operator determining whether the confidence
%   intervals will be calculated. Default value is FALSE.
if nargin < 4
    estimate_conf = false;
end
if nargin < 3
    write = false;
end
if nargin < 2
    figures = true;
end

% Load the result with minimal fval
results_dir = [ 'results/' analysis_name ];

display(results_dir)

files = dir([ results_dir '/INST_emu_*.mat' ]);
nresults = size(files, 1);
display(nresults)

fvalues = zeros(1, nresults);
xvalues = cell(1, nresults);
for ifile = 1:nresults
    load([results_dir '/' files(ifile).name], 'fval', 'x');
    fvalues(ifile) = fval;
    xvalues{ifile} = x;
end

% hist(fvalues);
[min_fval, imin_res] = min(fvalues);

%%

disp(files(imin_res));
load([results_dir '/' files(imin_res).name]);

x = xvalues{imin_res}; % best parameters x
if figures
    figure;
    plot_result(x, emumodel, ldata, ldata.src_isos);
end
fprintf('Objective of best solution: %g\n', min_fval);
fprintf('RMS standardized error of best solution: %g\n', ...
    sqrt(min_fval / sum(cellfun(@numel, ldata.mids))));
nz_isos = ldata.src_isos;
nrxn = size(emumodel.S, 2);
nmetab = length(emumodel.metabs);
nnz_frac = length(nz_isos);
nlmetab = length(ldata.metabs);

vfor = x(1:(nrxn - nnz(emumodel.reverse)));
vrev = zeros(nrxn - nnz(emumodel.reverse), 1);
vrev(emumodel.reverse) = x((nrxn - nnz(emumodel.reverse) + 1):nrxn);

opt.vnet = vfor - vrev;
opt.vxch = min(vfor, vrev);

opt.conc0 = x(nrxn + (1:nmetab));
opt.nz_fracs = x(nrxn + nmetab + (1:nnz_frac));

if exist('nlmetab', 'var')
    opt.dilutions = x(nrxn + nmetab + nnz_frac + (1:nlmetab));
end

opt.nparam = length(x);

if write
    write_results(analysis_name, opt, emumodel, ldata, nz_isos);
end

if estimate_conf
    % estimate confidence interval
    [~, result_name] = fileparts(files(imin_res).name);
    [~, opt_conf_int] = calc_conf_int(analysis_name, result_name, emumodel, ldata, nz_isos);
    
    if ~isempty(opt_conf_int)
        write_results_conf_int(analysis_name, opt, emumodel, ldata, nz_isos, opt_conf_int);
    end
end
