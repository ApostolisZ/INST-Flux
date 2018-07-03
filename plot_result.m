function plot_result(x, emumodel, ldata, nz_isos)
%PLOT_RESULT plots the MIDs and experimental data
%   PLOT_RESULT(X,EMUMODEL,LDATA,NZ_ISOS) generates a figure containing the
%   experimental data and the simulated data for the decision variables
%   contained in X
%   X is a vector containing the decision variables used in the model
%   EMUMODEL is a structure containing information on the EMU reactions in
%   the network
%   LDATA is a structure containing the labeling information on labeled
%   metabolites
%   NZ_ISOS contains information on source metabolites in the network

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

mids = simulate_emumodel(emumodel, v, conc, nz_isos, nz_fracs, src_metabs, tfinal);

nemu = size(ldata.emus, 1);
nplot = ceil(sqrt(nemu));
mplot = ceil(nemu / nplot);
MIDS = {};
MIDS{1,1} = 'Ion';
MIDS{1,2} = 'EMU';
MIDS{1,3} = 'Time (min)';
for i = 1:numel(ldata.t)
    MIDS{1,i+3} = ldata.t(i);
end
fin = 1;
for iemu = 1:nemu
    
    subplot(mplot, nplot, iemu)

    sim_mid = calc_sim_mid(mids, iemu, emumodel, ldata, conc, dilutions);

    plot(ldata.t, sim_mid');
    
    hold on
    set(gca, 'ColorOrderIndex', 1);
    if size(ldata.raw, 2) > 1
        errorbar(repmat(ldata.t, size(ldata.mids{iemu}, 1), 1)', ldata.mids{iemu}', ldata.mids_std{iemu}', 'x');
    else
        plot(ldata.t, ldata.mids{iemu}', 'x')
    end
    title(strrep(ldata.emus(iemu), '_', '\_'))
    xlabel('Time (min)')
    ylabel('Fractional abundance')
    
    xlim([min(ldata.t) max(ldata.t)]);
    yl = ylim;
    ylim([max(0, yl(1)) min(1, yl(2))]);
    
    nmids = size(sim_mid, 1) - 1;
    Mstr = [repmat(['M'], nmids + 1,1) num2str((0:nmids)')];
%     axisLimits = axis;
%     axis([0 15 axisLimits(3:4)]);
    legend([repmat(['M + '], nmids + 1,1) num2str((0:nmids)')])
    hold off
    start = fin+1;
    fin = start+nmids;
    MIDS{start,1} = ldata.ions{iemu};
    MIDS{start,2} = ldata.emus{iemu};
    rand_mid = sim_mid; %+ rand*0.00125.*sim_mid - rand*0.00125.*sim_mid;
    for i = start:fin
    MIDS{i,3} = Mstr(i-start+1,:);
        for j = 1:size(sim_mid,2)
            MIDS{i,j+3} =10000*rand*rand_mid(i-start+1,j);
        end
    end
%     start = fin+1;
%     fin = start+nmids;
%     rand_mid = sim_mid + rand*0.0125.*sim_mid - rand*0.0125.*sim_mid;
%     for i = start:fin
%     MIDS{i,3} = Mstr(i-start+1,:);
%         for j = 1:size(sim_mid,2)
%             MIDS{i,j+3} = 10000*rand*rand_mid(i-start+1,j);
%         end
%     end
end

% for i = 1:8
%     h = figure(i);
%     saveas(h, ['figure' num2str(h)], 'png');
% end

end