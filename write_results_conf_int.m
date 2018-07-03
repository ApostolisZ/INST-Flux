function write_results_conf_int(analysis_name, opt, emumodel, ldata, nz_isos, opt_conf_int)

fid = fopen([ analysis_name '/results-' analysis_name '-interval.csv' ], 'w');

fprintf(fid, 'Fluxes\n');
fprintf(fid, 'Reaction,Net flux,Net flux lb,Net flux ub,Exchange flux,Exchange flux lb Exchange flux ub\n');
for irxn = 1:length(opt.vnet)
    rxnstring = print_rxn(emumodel.S(:, irxn), emumodel.metabs, emumodel.reverse(irxn));
    fprintf(fid, '%s,%g,%g,%g,%g,%g,%g\n', rxnstring, opt.vnet(irxn), opt_conf_int.vnet_lb(irxn), opt_conf_int.vnet_ub(irxn), opt.vxch(irxn), opt_conf_int.vxch_lb(irxn), opt_conf_int.vxch_ub(irxn));
end

fprintf(fid, '\nMetabolite Concentrations\n');
fprintf(fid, 'Metabolite,Concentration,Concentration lb,Concentration ub\n');
for imetab = 1:length(opt.conc0)
    fprintf(fid, '%s,%g,%g,%g\n', emumodel.metabs{imetab}, opt.conc0(imetab), opt_conf_int.conc0_lb(imetab), opt_conf_int.conc0_ub(imetab));
end

fprintf(fid, '\nIsotopomer Fractions\n');
fprintf(fid, 'Isotopomer,Fraction,Fraction lb,Fraction ub\n');
for inz_iso = 1:length(opt.nz_fracs)
    fprintf(fid, '%s,%g,%g,%g\n', nz_isos{inz_iso}, opt.nz_fracs(inz_iso), opt_conf_int.nz_fracs_lb(inz_iso), opt_conf_int.nz_fracs_ub(inz_iso));
end

fprintf(fid, '\nDilution Factors\n');
fprintf(fid, 'Metabolite,Dilution factor,Dilution factor lb,Dilution factor ub\n');
for idilution = 1:length(opt.dilutions)
    fprintf(fid, '%s,%g,%g,%g\n', ldata.metabs{idilution}, opt.dilutions(idilution), opt_conf_int.dilutions_lb(idilution), opt_conf_int.dilutions_ub(idilution));
end

fclose(fid);

end

function rxnstring = print_rxn(Scol, metabs, reverse)
rxnstring = '';
if any(Scol < 0)
    for imetab = reshape(find(Scol < 0), 1, [])
        rxnstring = [ rxnstring sprintf('%g %s + ', -full(Scol(imetab)), metabs{imetab}) ];
    end
    rxnstring = rxnstring(1:(end - 2));
end

if reverse
    rxnstring = [ rxnstring '= ' ];
else
    rxnstring = [ rxnstring '-> ' ];
end

if any(Scol > 0)
    for imetab = reshape(find(Scol > 0), 1, [])
        rxnstring = [ rxnstring sprintf('%g %s + ', full(Scol(imetab)), metabs{imetab}) ];
    end
    rxnstring = rxnstring(1:(end - 2));
end
end