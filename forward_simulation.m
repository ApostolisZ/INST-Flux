function sim_mids = forward_simulation(analysis_name,x)
NETWORK_FNAME = 'network.csv';
ATOM_SEP = '_';

if nargin < 2
    init_fracs_zero = true;
end

TOL = 1e-10;
MIN_VAL = 1e-6;
MAX_VAL = 1e6;
MAX_V0_MULT = 2;

disp(analysis_name);

ldata = read_ldata_forward_sim(analysis_name);

if isempty(ldata)
    return
end
nlmetab = length(ldata.metabs);
bnd_metabs = unique(cellfun(@(x) x(1:(strfind(x, ATOM_SEP) - 1)), ldata.src_isos, 'UniformOutput', false));
emumodel = generate_emumodel([analysis_name '/' NETWORK_FNAME], ldata.emus, bnd_metabs);  

[emus_blocks_sorted,emus_blocks_order] = sort(emumodel.emu_block);
new_EMU_list = emumodel.emus(emus_blocks_order);
emus_blocks_sorted = emus_blocks_sorted';
Emus = {new_EMU_list,emus_blocks_sorted};



% Determine set of "non-zero" isotopomers, i.e. isotopomers whose initial
% fractional composition of the pool may be non-zero.
nz_isos = ldata.src_isos;

% Read pool concentrations
cdata = read_conc_forward_sim(analysis_name, emumodel, bnd_metabs);  
conc = cdata.x;

% read fluxes
v = emumodel.bfixed;
% read substrate
nz_isos = ldata.src_isos;
nz_fracs = ldata.src_fracs;
n = regexp(ldata.src_isos, '(?<metab>\S+)_(?<iso>\S+)', 'names');
src_metabs = unique(cellfun(@(x) x.metab, n, 'UniformOutput', false));
tfinal = ldata.t(end);
sim_mids = plot_result_forward_sim(x, emumodel, ldata, nz_isos);

end



