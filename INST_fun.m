function y = INST_fun(x, emumodel, ldata)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

% EPS = 1e-6;

global mids;

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
    
    sim_mids = cell(length(emu_inds), 1);
    for i = 1:length(ldata.emus)
        emu = ldata.emus{i};
        k = strfind(emu, ATOM_SEP);
        metab = emu(1:(k(end) - 1));
        atom_nums = emu((k(end) + 1):end);
        dilution = dilutions(i);
        
        mid_inds = find(~cellfun('isempty', regexp(emumodel.emus, [ '^' metab '(\' COMPART_SEP '\w+|)\' ATOM_SEP atom_nums '$' ])));
        sim_mid = zeros(size(ldata.mids{i}));
        conc = zeros(length(mid_inds), 1);
        for imid = reshape(mid_inds, 1, [])
            k = strfind(emumodel.emus{imid}, ATOM_SEP);
            compart_metab = emumodel.emus{imid}(1:(k(end) - 1));
            conc(imid) = conc(strcmp(compart_metab, emumodel.metabs));
            sim_mid(2:end, :) = sim_mid(2:end, :) ...
                + (deval(mids{imid}{1}, ldata.t, mids{imid}{2}{1})) * conc(imid);
        end
        sim_mid(2:end, :) = sim_mid(2:end, :) / sum(conc) * dilution;
        sim_mid(1, :) = 1 - sum(sim_mid(2:end, :), 1);
        
        sim_mids{i} = sim_mid;
    end
    
    y = reshape(cell2mat(sim_mids), [], 1);
catch err
    disp(getReport(err));
    save;
    
    y = nan(sum(cellfun(@numel, ldata.mids)), 1);
end
end

