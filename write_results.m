function write_results(analysis_name, opt, emumodel, ldata, nz_isos)
%WRITE_RESULTS generates an output file 
%   WRITE_RESULTS(ANALYSIS_NAME,OPT,EMUMODEL,LDATA,NZ_ISOS) generates an
%   output file "results-ANALYSIS_NAME.csv" in the directory ANALYSIS_NAME, containing
%   the information in the result structure OPT. This describes the EMU
%   model EMUMODEL with labeling data LDATA and source metabolites
%   NZ_ISOS. 
%
%   ANALYSIS_NAME, OPT are strings containing directory and file names,
%   EMUMODEL is a model structure
%   LDATA is a structure containing the labeling data
%   NZ_ISOS contains information on source metabolites in the network

fid = fopen([ analysis_name '/results-' analysis_name '.csv' ], 'w');

fprintf(fid, 'Fluxes\n');
fprintf(fid, 'Reaction,Net flux,Exchange flux\n');
for irxn = 1:length(opt.vnet)
    rxnstring = print_rxn(emumodel.S(:, irxn), emumodel.metabs, emumodel.reverse(irxn));
    fprintf(fid, '%s,%g,%g\n', rxnstring, opt.vnet(irxn), opt.vxch(irxn));
end

fprintf(fid, '\nMetabolite Concentrations\n');
fprintf(fid, 'Metabolite,Concentration\n');
for imetab = 1:length(opt.conc0)
    fprintf(fid, '%s,%g\n', emumodel.metabs{imetab}, opt.conc0(imetab));
end

fprintf(fid, '\nIsotopomer Fractions\n');
fprintf(fid, 'Isotopomer,Fraction\n');
for inz_iso = 1:length(opt.nz_fracs)
    fprintf(fid, '%s,%g\n', nz_isos{inz_iso}, opt.nz_fracs(inz_iso));
end

fprintf(fid, '\nDilution Factors\n');
fprintf(fid, 'Metabolite,Dilution factor\n');
for idilution = 1:length(opt.dilutions)
    fprintf(fid, '%s,%g\n', ldata.metabs{idilution}, opt.dilutions(idilution));
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