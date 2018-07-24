function conc = read_conc_forward_sim(analysis_dir, emumodel, bnd_metabs)
%READ_CDATA read the laeling data and correct for natural abundance of
%stable isotopes.
%   CONC = READ_CDATA(ANALYSIS_DIR, EMUMODEL, BND_METABS) navigates 
%   to ANALYSIS_DIR and reads the concentration data in the data file
%   conc.csv. 
%
%   CONC : a structure including information about the measured pool sizes
%   of network metabolites with the following fields
%   inds : the indeces of the measured metabolites in the metabolites list
%   metabs : the metabolites with measured pool sizes
%   x : the metabolite concentrations
%   x_std : the metabolite concentration standard error
%   ANALYSIS_DIR is the directory including the data files.
%   CORRECT is a boolean variable, indicating whenter the eperimental data
%   will be corrected based on natural abundance.
DATA_FNAME = 'pool_forward.csv';
COMPART_SEP = '.';

fid = fopen([ analysis_dir '/' DATA_FNAME ]);
if fid == -1
    conc.inds = zeros(0, 1);    
    conc.x = zeros(0, 1);
    conc.x_std = zeros(0, 1);
else
    C = textscan(fid, '%s%f%f', 'Delimiter', ',', 'Headerlines', 1);
    fclose(fid);
    
    metabs_without_compart = cellfun(@(x) x(1), ...
        cellfun(@(x) strsplit(x, COMPART_SEP), emumodel.metabs, 'UniformOutput', false))';
    conc.inds = cellfun(@(x) find(strcmp(x, metabs_without_compart)), C{1}, ...
        'UniformOutput', false);
    conc.metabs = C{1};
    conc.x = C{2};
        
    conc.x_std = nan(size(conc.x));
    conc.x_std(1:length(C{3})) = C{3};
    conc.x_std(isnan(conc.x_std)) = 0;
end

for i = 1:length(bnd_metabs)
    imetab = find(strcmp(bnd_metabs{i}, emumodel.metabs));
    if isempty(conc.inds) || ~any(cellfun(@(x) ismember(imetab, x), conc.inds))
        conc.metabs{end + 1} = bnd_metabs{i};
        conc.inds{end + 1} = imetab;
        conc.x(end + 1) = 1;
        conc.x_std(end + 1) = 0;
    end
end

end