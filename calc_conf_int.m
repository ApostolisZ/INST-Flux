function [xvalues, conf_int] = calc_conf_int(analysis_name, result_name, emumodel, ldata, nz_isos)
% [XVALUES, CONF_INT] = CALC_CONF_INT(ANALYSIS_NAME, RESULT_NAME, EMUMODEL,
% LDATA, NZ_ISOS) calculates the confidence intervals around an optimal
% solution, by using the Monte Carlo method. The optimal solutions with
% added noise need to be calculated previously, by utilizing the function
% calc_ci
% XVALUES is the parameter values calculated based on the noisy data
% CONF_INT is a structure which contains the upper and lower bound for 95%
% confidence.
% ANALYSIS_NAME is the name of the folder containing the data set
% RESULT_NAME is the name of the calculated optimal result
% EMUMODEL is a structure containing the emu network
% LDATA is a structure containing information on the labeling data
% NZ_ISOS is a celll array containing the labeled isotopomers

results_dir = [ 'results/' analysis_name '/' result_name ];

files = dir([ results_dir '/ci_*.mat' ]);
nresults = size(files, 1);
display(nresults)

xvalues = cell(1, nresults);
for ifile = 1:nresults
    load([results_dir '/' files(ifile).name], 'fval', 'x');
    xvalues{ifile} = x;
end
xvalues = cell2mat(xvalues);

nrxn = size(emumodel.S, 2);
nmetab = length(emumodel.metabs);
nnz_frac = length(nz_isos);
nlmetab = length(ldata.metabs);

%
% Estimate confidence interval
%
p = [0.025, 0.975];

vfor = xvalues(1:(nrxn - nnz(emumodel.reverse)), :);
vrev = zeros(nrxn - nnz(emumodel.reverse), size(xvalues, 2));
vrev(emumodel.reverse, :) = xvalues((nrxn - nnz(emumodel.reverse) + 1):nrxn, :);

vnet = vfor - vrev;
vxch = min(vfor, vrev);

conc0 = xvalues(nrxn + (1:nmetab), :);
nz_fracs = xvalues(nrxn + nmetab + (1:nnz_frac), :);

dilutions = xvalues(nrxn + nmetab + nnz_frac + (1:nlmetab), :);

conf_int_vnet = quantile(vnet, p, 2);
conf_int_vxch = quantile(vxch, p, 2);
conf_int_conc0 = quantile(conc0, p, 2);
conf_int_nz_fracs = quantile(nz_fracs, p, 2);
conf_int_dilutions = quantile(dilutions, p, 2);

conf_int.vnet_lb = conf_int_vnet(:, 1);
conf_int.vnet_ub = conf_int_vnet(:, 2);
conf_int.vxch_lb = conf_int_vxch(:, 1);
conf_int.vxch_ub = conf_int_vxch(:, 2);
conf_int.conc0_lb = conf_int_conc0(:, 1);
conf_int.conc0_ub = conf_int_conc0(:, 2);
conf_int.nz_fracs_lb = conf_int_nz_fracs(:, 1);
conf_int.nz_fracs_ub = conf_int_nz_fracs(:, 2);
conf_int.dilutions_lb = conf_int_dilutions(:, 1);
conf_int.dilutions_ub = conf_int_dilutions(:, 2);
end