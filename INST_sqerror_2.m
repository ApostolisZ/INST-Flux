function [sqerror, dsqerror] = INST_sqerror_2(x, emumodel, ldata, nz_isos)
%INST_SQERROR calculates the value of the objective function
%   [SQERROR, DSQERROR] = INST_SQERROR(X,EMUMODE,LDATA,NZ_ISOS)
%   calculates the value of the objective function used in INST-FLUX. The
%   objective funtion is as defined on Young et al. (Biotech. and Bioeng., 99(3):686-99,
%   2007). 
%   SQERROR is the value of the objective function
%   DSQERROR is the jacobian of the objective function
%	X is a vector of the decision variables used to calculate the objective
%	function
%   EMUMODEL is a structure containing information on the EMU reactions in
%   the network
%   LDATA is a structure containing the labeling information on labeled
%   metabolites
%   NZ_ISOS contains information on source metabolites in the network
%   MIDS is the solution structure produced when solving the ODEs
%   containing information on the labeling of network EMUs

nmetab = size(emumodel.S, 1);
nrxn = size(emumodel.S, 2);
nnz_iso = length(nz_isos);
nlmetab = length(ldata.metabs);

x = reshape(x, [], 1);
v = x(1:nrxn);
conc = x(nrxn + (1:nmetab));
nz_fracs = x(nrxn + nmetab + (1:nnz_iso));
dilutions = x(nrxn + nmetab + nnz_iso + (1:nlmetab));
tfinal = ldata.t(end);

n = regexp(ldata.src_isos, '(?<metab>\S+)_(?<iso>\S+)', 'names');
src_metabs = unique(cellfun(@(x) x.metab, n, 'UniformOutput', false));

try
    mids = simulate_emumodel(emumodel, v, conc, nz_isos, nz_fracs, src_metabs, tfinal);
    
    sqerror = 0.0;
    for i = 1:length(ldata.emus)
        sim_mid = calc_sim_mid(mids, i, emumodel, ldata, conc, dilutions);
        meas_mid = ldata.mids{i};
        meas_mid_std = ldata.mids_std{i};
        
        sim_mid = sim_mid(~isnan(meas_mid));
        meas_mid_std = meas_mid_std(~isnan(meas_mid));
        meas_mid = meas_mid(~isnan(meas_mid));
        
        sqerror = sqerror + ...
            sum(((sim_mid - meas_mid) ./ meas_mid_std).^2);
    end
    
    if nargout > 1
        dsqerror = INST_sqerror_grad(x, emumodel, ldata, nz_isos, mids);
    end
catch err
    disp(getReport(err));
    
    sqerror = nan;
    dsqerror = nan(nmetab + nrxn + nnz_iso + nlmetab, 1);
end

end

