function [sim_mid, metab, dilution, mid_inds, sim_mids] = calc_sim_mid(mids, iemu, emumodel, ldata, conc0, dilutions)
%[SIM_MID, METAB, DILUTION, MIN_INDS, SIM_MIDS] = CALC_SIM_MID(MIDS, IEMU, EMUMODEL, LDATA, CONC0, DILUTIONS)
% calculates the MIDs for emu IEMU, based on the ODE solution structure
% MIDS. 
%   MIDS is the solution structure produced when solving the ODEs
%   containing information on the labeling of network EMUs
%   IEMU is the emu index for which the MIDs are calculated
%   LDATA is a structure containing the labeling data
%   CONC0 is a vector containing the metabolite concentrations
%   DILUTIONS is a vector that contains the dilution factors for the
%   labeled metabolites.


ATOM_SEP = '_';
COMPART_SEP = '.';

emu = ldata.emus{iemu};
k = strfind(emu, ATOM_SEP);
metab = emu(1:(k(end) - 1));
atom_nums = emu((k(end) + 1):end);
dilution = dilutions(strcmp(metab, ldata.metabs));

mid_inds = find(~cellfun('isempty', regexp(emumodel.emus, [ '^' metab '(\' COMPART_SEP '\w+|)\' ATOM_SEP atom_nums '$' ])));
metabs = cell(length(mid_inds), 1);
conc = zeros(length(mid_inds), 1);

sim_mids = zeros(length(mid_inds), size(ldata.mids{iemu}, 1), size(ldata.mids{iemu}, 2));
for imid = 1:length(mid_inds)
    mid_ind = mid_inds(imid);
    k = strfind(emumodel.emus{mid_ind}, ATOM_SEP);
    metabs{imid} = emumodel.emus{mid_ind}(1:(k(end) - 1));
    conc(imid) = conc0(strcmp(metabs{imid}, emumodel.metabs));
    
    sim_mids(imid, 2:end, :) = ...
        + (deval(mids{mid_ind}{1}, ldata.t, mids{mid_ind}{2}{1})) * dilution;
    sim_mids(imid, 1, :) = 1 - sum(sim_mids(imid, 2:end, :), 2);
end

sim_mid = squeeze(sum(sim_mids .* repmat(conc, 1, size(sim_mids, 2), size(sim_mids, 3)), 1)) ...
    ./ sum(conc);
end