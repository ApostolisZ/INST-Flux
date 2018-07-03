function grad = INST_grad(x, emumodel, ldata)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

global mids;
global dmids;

nmetab = size(emumodel.S, 1);
nrxn = size(emumodel.S, 2);
nsrc_iso = length(ldata.src_isos);
nlmetab = length(ldata.metabs);

x = reshape(x, [], 1);
v = x(1:nrxn);
conc = x(nrxn + (1:nmetab));
src_fracs = x(nrxn + nmetab + (1:nsrc_iso));
dilutions = x(nrxn + nmetab + nsrc_iso + (1:nlmetab));
tfinal = ldata.t(end);

try
    simulate_emumodel(emumodel, v, conc, ldata.src_isos, src_fracs, tfinal, 0); 
    
    dsim_mids = cell(length(emu_inds), nrxn + nmetab + nsrc_iso + nlmetab);
    for i = 1:length(ldata.emus)
        emu = ldata.emus{i};
        k = strfind(emu, ATOM_SEP);
        metab = emu(1:(k(end) - 1));
        atom_nums = emu((k(end) + 1):end);
        dilution = dilutions(i);
        
        mid_inds = find(~cellfun('isempty', regexp(emumodel.emus, [ '^' metab '(\' COMPART_SEP '\w+|)\' ATOM_SEP atom_nums '$' ])));
        sim_mids = zeros(length(mid_inds), size(ldata.mids{i}, 1), size(ldata.mids{i}, 2));
        metabs = cell(length(mid_inds), 1);
        mids_concs = zeros(length(mid_inds), 1);
        for imid = 1:length(mid_inds)
            mid_ind = mid_inds(imid);
            k = strfind(emumodel.emus{mid_ind}, ATOM_SEP);
            metabs{imid} = emumodel.emus{mid_ind}(1:(k(end) - 1));
            mids_concs(imid) = conc(strcmp(metabs{imid}, emumodel.metabs));
            sim_mids(imid, 2:end, :) = ...
                + (deval(mids{mid_ind}{1}, ldata.t, mids{mid_ind}{2}{1}));
            sim_mids(imid, 1, :) = 1 - sum(sim_mids(imid, 2:end, :), 2);
        end
        sim_mid = squeeze(sum(sim_mids .* repmat(mids_concs, 1, size(sim_mids, 2), size(sim_mids, 3)), 1)) ...
            ./ sum(mids_concs);
        
        for irxn = 1:nrxn
            simulate_emumodel(emumodel, v, conc, ldata.src_isos, src_fracs, tfinal, irxn);
            
            dsim_mid = zeros(size(ldata.mids{i}));
            for imid = 1:length(mid_inds)
                mid_ind = mid_inds(imid);
                dsim_mid(2:end, :) = dsim_mid(2:end, :) ...
                    + (deval(dmids{mid_ind}{1}, ldata.t, dmids{mid_ind}{2}{1})) * mids_concs(imid);
            end
            dsim_mid(2:end, :) = dsim_mid(2:end, :) / sum(mids_concs) * dilution;
            dsim_mid(1, :) = -sum(dsim_mid(2:end, :), 1);
            
            dsim_mids{i, irxn} = dsim_mid;
        end
        
        for imetab = 1:nmetab
            simulate_emumodel(emumodel, v, conc, ldata.src_isos, src_fracs, tfinal, nrxn + imetab);
            
            dsim_mid = zeros(size(ldata.mids{i}));
            for imid = 1:length(mid_inds)
                mid_ind = mid_inds(imid);
                dsim_mid(2:end, :) = dsim_mid(2:end, :) ...
                    + (deval(dmids{mid_ind}{1}, ldata.t, dmids{mid_ind}{2}{1})) * mids_concs(imid);
            end
            dsim_mid(2:end, :) = dsim_mid(2:end, :) / sum(mids_concs) * dilution;
            dsim_mid(1, :) = -sum(dsim_mid(2:end, :), 1);
            
            if any(strcmp(emumodel.metabs(imetab), metabs))
                dsim_mid = dsim_mid ...
                    + (squeeze(sim_mids(strcmp(emumodel.metabs(7), metabs), :, :)) - sim_mid) / sum(mids_concs);
            end
            
            dsim_mids{i, nrxn + imetab} = dsim_mid;
        end
        
        for isrc_iso = 1:nsrc_iso
            simulate_emumodel(emumodel, v, conc, ldata.src_isos, src_fracs, tfinal, nrxn + nmetab + isrc_iso);
            
            dsim_mid = zeros(size(ldata.mids{i}));
            for imid = 1:length(mid_inds)
                mid_ind = mid_inds(imid);
                dsim_mid(2:end, :) = dsim_mid(2:end, :) ...
                    + (deval(dmids{mid_ind}{1}, ldata.t, dmids{mid_ind}{2}{1})) * mids_concs(imid);
            end
            dsim_mid(2:end, :) = dsim_mid(2:end, :) / sum(mids_concs) * dilution;
            dsim_mid(1, :) = -sum(dsim_mid(2:end, :), 1);
            
            dsim_mids{i, nrxn + nmetab + isrc_iso} = dsim_mid;
        end
        
        idilution = find(strcmp(metab, ldata.metabs));
        dsim_mid = zeros(size(ldata.mids{i}));
        dsim_mid(2:end, :) = sim_mid(2:end, :) / dilution;
        dsim_mid(1, :) = -sum(dsim_mid, 1);
        
        dsim_mids{i, nrxn + nmetab + nsrc_iso + idilution} = dsim_mid;
        for jdilution = setdiff(1:nlmetab, idilution)
            dsim_mids{i, nrxn + nmetab + nsrc_iso + jdilution} = zeros(size(dsim_mid));
        end
    end
    
    grad = reshape(cell2mat(dsim_mids), [], nrxn + nmetab + nsrc_iso + nlmetab);
catch err
    disp(getReport(err));
    
    grad = nan(sum(cellfun(@numel, ldata.mids)), nmetab + nrxn + nsrc_iso + nlmetab);
end
end

