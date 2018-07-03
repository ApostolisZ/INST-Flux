function d2sqerror = INST_sqerror_2grad(x,lambda, emumodel, ldata, nz_isos, mids)
% function [dsqerror,d2sqerror] = INST_sqerror_2grad(x, emumodel, ldata, nz_isos, mids)
%INST_SQERROR_GRAD Summary of this function goes here
%   Detailed explanation goes here
global Jac
ATOM_SEP = '_';


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


Sigma = [];

if nargin < 4
    mids = [];
end



try
    if isempty(mids)
        mids = simulate_emumodel(emumodel, v, conc, nz_isos, nz_fracs, src_metabs, tfinal); 
    end
        F = [];
    dsqerror = zeros(nmetab + nrxn + nnz_iso + nlmetab, 1);       
    for irxn = 1:nrxn

        dmids = simulate_demumodel(emumodel, v, conc, nz_isos, nz_fracs, src_metabs, tfinal, mids, irxn);        
        
        for i = 1:length(ldata.emus)
            [sim_mid, ~, dilution, mid_inds] = calc_sim_mid(mids, i, emumodel, ldata, conc, dilutions);
            
            metabs = cell(length(mid_inds), 1);
            mids_concs = zeros(length(mid_inds), 1);            
            dsim_mid = zeros(size(ldata.mids{i}));
            for imid = 1:length(mid_inds)
                mid_ind = mid_inds(imid);
                k = strfind(emumodel.emus{mid_ind}, ATOM_SEP);
                metabs{imid} = emumodel.emus{mid_ind}(1:(k(end) - 1));
                mids_concs(imid) = conc(strcmp(metabs{imid}, emumodel.metabs));
                
                dsim_mid(2:end, :) = dsim_mid(2:end, :) ...
                    + (deval(dmids{mid_ind}{1}, ldata.t, dmids{mid_ind}{2}{1})) * mids_concs(imid);
            end
                                
            dsim_mid(2:end, :) = dsim_mid(2:end, :) / sum(mids_concs) * dilution;
            dsim_mid(1, :) = -sum(dsim_mid(2:end, :), 1);
            
            meas_mid = ldata.mids{i};
            meas_mid_std = ldata.mids_std{i};
            
            sim_mid = sim_mid(~isnan(meas_mid));
            dsim_mid = dsim_mid(~isnan(meas_mid));
            meas_mid_std = meas_mid_std(~isnan(meas_mid));
            meas_mid = meas_mid(~isnan(meas_mid));
            F = [F;sim_mid - meas_mid];
            J{i,irxn} = dsim_mid;
            SIgma = [Sigma;meas_mid_std];
%             H(i,:) = dsim_mid;
            dsqerror(irxn) = dsqerror(irxn) + ...
                2 * sum(sum((sim_mid - meas_mid) ./ (meas_mid_std.^2) .* dsim_mid));

        end
    end
        
    for imetab = 1:nmetab
        dmids = simulate_demumodel(emumodel, v, conc,  nz_isos, nz_fracs, src_metabs, tfinal, mids, nrxn + imetab);
        for i = 1:length(ldata.emus)
            [sim_mid, ~, dilution, mid_inds, sim_mids] = calc_sim_mid(mids, i, emumodel, ldata, conc, dilutions);
            
            metabs = cell(length(mid_inds), 1);
            mids_concs = zeros(length(mid_inds), 1);                        
            dsim_mid = zeros(size(ldata.mids{i}));
            for imid = 1:length(mid_inds)
                mid_ind = mid_inds(imid);
                k = strfind(emumodel.emus{mid_ind}, ATOM_SEP);
                metabs{imid} = emumodel.emus{mid_ind}(1:(k(end) - 1));
                mids_concs(imid) = conc(strcmp(metabs{imid}, emumodel.metabs));
                
                dsim_mid(2:end, :) = dsim_mid(2:end, :) ...
                    + (deval(dmids{mid_ind}{1}, ldata.t, dmids{mid_ind}{2}{1})) * mids_concs(imid);
            end
            
            dsim_mid(2:end, :) = dsim_mid(2:end, :) / sum(mids_concs) * dilution;
            dsim_mid(1, :) = -sum(dsim_mid(2:end, :), 1);
            
            if any(strcmp(emumodel.metabs(imetab), metabs))
                dsim_mid = dsim_mid ...
                    + (squeeze(sim_mids(strcmp(emumodel.metabs(imetab), metabs), :, :)) - sim_mid) / sum(mids_concs);
            end
            
            meas_mid = ldata.mids{i};
            meas_mid_std = ldata.mids_std{i};            
            
            sim_mid = sim_mid(~isnan(meas_mid));
            dsim_mid = dsim_mid(~isnan(meas_mid));
            meas_mid_std = meas_mid_std(~isnan(meas_mid));
            meas_mid = meas_mid(~isnan(meas_mid));
            J{i,irxn+imetab} = dsim_mid;
%             SIgma = [Sigma;meas_mid_std];
            dsqerror(nrxn + imetab) = dsqerror(nrxn + imetab) + ...
                2 * sum(sum((sim_mid - meas_mid) ./ (meas_mid_std.^2) .* dsim_mid));

        end
    end
        
    for inz_iso = 1:nnz_iso
        dmids = simulate_demumodel(emumodel, v, conc,  nz_isos, nz_fracs, src_metabs, tfinal, mids, nrxn + nmetab + inz_iso);
        for i = 1:length(ldata.emus)
            [sim_mid, ~, dilution, mid_inds] = calc_sim_mid(mids, i, emumodel, ldata, conc, dilutions);
                        
            metabs = cell(length(mid_inds), 1);
            mids_concs = zeros(length(mid_inds), 1);            
            dsim_mid = zeros(size(ldata.mids{i}));
            for imid = 1:length(mid_inds)
                mid_ind = mid_inds(imid);
                k = strfind(emumodel.emus{mid_ind}, ATOM_SEP);
                metabs{imid} = emumodel.emus{mid_ind}(1:(k(end) - 1));
                mids_concs(imid) = conc(strcmp(metabs{imid}, emumodel.metabs));
                
                dsim_mid(2:end, :) = dsim_mid(2:end, :) ...
                    + (deval(dmids{mid_ind}{1}, ldata.t, dmids{mid_ind}{2}{1})) * mids_concs(imid);
            end
                       
            dsim_mid(2:end, :) = dsim_mid(2:end, :) / sum(mids_concs) * dilution;
            dsim_mid(1, :) = -sum(dsim_mid(2:end, :), 1);
            
            meas_mid = ldata.mids{i};
            meas_mid_std = ldata.mids_std{i};
            
            sim_mid = sim_mid(~isnan(meas_mid));
            dsim_mid = dsim_mid(~isnan(meas_mid));
            meas_mid_std = meas_mid_std(~isnan(meas_mid));
            meas_mid = meas_mid(~isnan(meas_mid));            
            J{i,irxn+imetab+inz_iso} = dsim_mid;
%             SIgma = [Sigma;meas_mid_std];
            dsqerror(nrxn + nmetab + inz_iso) = dsqerror(nrxn + nmetab + inz_iso) + ...
                2 * sum(sum((sim_mid - meas_mid) ./ (meas_mid_std.^2) .* dsim_mid));

        end
    end
        
    for i = 1:length(ldata.emus)
        [sim_mid, metab, dilution] = calc_sim_mid(mids, i, emumodel, ldata, conc, dilutions);        
        idilution = find(strcmp(metab, ldata.metabs));
        dsim_mid = zeros(size(ldata.mids{i}));
        dsim_mid(2:end, :) = sim_mid(2:end, :) / dilution;
        dsim_mid(1, :) = -sum(dsim_mid, 1);
        
        meas_mid = ldata.mids{i};
        meas_mid_std = ldata.mids_std{i};
        
        sim_mid = sim_mid(~isnan(meas_mid));
        dsim_mid = dsim_mid(~isnan(meas_mid));
        meas_mid_std = meas_mid_std(~isnan(meas_mid));
        meas_mid = meas_mid(~isnan(meas_mid));
        J{i,irxn+imetab+nnz_iso+i} = dsim_mid;
%         SIgma = [Sigma;meas_mid_std];
        dsqerror(nrxn + nmetab + nnz_iso + idilution) = dsqerror(nrxn + nmetab + nnz_iso + idilution) + ...
                2 * sum(sum((sim_mid - meas_mid) ./ (meas_mid_std.^2) .* dsim_mid));

        
    end
    Jm = cell2mat(J);
    Jac = Jm;
     d2sqerror = Jm'*diag(SIgma)*Jm;
catch err
    disp(getReport(err));
    
    dsqerror = nan(nmetab + nrxn + nnz_iso + nlmetab, 1);
end

end

