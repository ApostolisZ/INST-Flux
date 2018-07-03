function data = read_ldata(analysis_dir, correct)
%READ_LDATA read the laeling data and correct for natural abundance of
%stable isotopes.
%   DATA = READ_LDATA(ANALYSIS_DIR, CORRECT) navigates to ANALYSIS_DIR and reads
%   the labeling state and fractions of the network substrate, as well as
%   the time-dependent labeling data of network metabolites. Metabolite labeling data are 
%   corrected for natural abundance of isotopic labels.
%   
%   DATA : a structure including labeling data of the substrate and
%   network metabolites with the following fields
%   t : the time points that measurements were taken
%   ions : the metabolite ions measured 
%   emus : the metabolite EMUs measured
%   raw : raw data for metabolite labeling
%   mids : corrected data for metabolites measured
%   mids_std : standard error for metabolites measured
%   metabs : measured metabolites
%   src_isos : source isotopes of the network
%   src_fracs : fractions of the total pool of the source isotopes
%   
%   ANALYSIS_DIR : the directory including the data files.
%   CORRECT is a boolean variable, indicating whenter the eperimental data
%   will be corrected based on natural abundance.

DATA_FNAME = 'labeling.csv';
SUBSTRATE_FNAME = 'substrate.csv';
ABUNDANCE_FNAME = 'abundance.csv';
ELE_COMP_FNAME  = 'element_comp.csv';
MIN_STD = 1e-3;     % minimum std for data point

ATOM_SEP = '_';

if nargin < 2
    correct = true;
end

[elements, abundances] = load_abundances(ABUNDANCE_FNAME);
ion_comps = load_compositions([analysis_dir '/' ELE_COMP_FNAME]);

A = importdata([analysis_dir '/' DATA_FNAME], ',');

data.t = A.data(1, :);          % t = time point array
head = A.textdata(2:end, 1);
data.ions = head(~cellfun('isempty', head));
head = A.textdata(2:end, 2);
data.emus = head(~cellfun('isempty', head));

data.raw = cell(length(data.emus), 1);
data.mids = cell(length(data.emus), 1);
data.mids_std = cell(length(data.emus), 1);
data.metabs = {};
for iion = 1:length(data.ions)
    ion_startrow = find(strcmp(data.ions(iion), A.textdata(:, 1)));
    if iion < length(data.emus)
        ion_endrow = find(strcmp(data.ions(iion + 1), A.textdata(:, 1))) - 1;
    else
        ion_endrow = length(A.textdata(:, 1));
    end
    
    rep_startrows = [ion_startrow + find(strcmp(A.textdata(ion_startrow:ion_endrow, 3)', 'M0')) - 1 ion_endrow + 1];
    for irep = 1:(length(rep_startrows) - 1)
        data.raw{iion, irep} = A.data(rep_startrows(irep):(rep_startrows(irep + 1) - 1), :);
    end
    
    k = strfind(data.emus{iion}, ATOM_SEP);
    data.metabs = union(data.emus{iion}(1:(k(end) - 1)), data.metabs);
end

% Correct for atmospheric abundances
corr_matrices = cell(size(data.ions));
for iion = 1:length(data.ions)
    ion_ind = find(strcmp(ion_comps.ions, data.ions{iion})); %just in case different order
    if isempty(ion_ind)
        error('Ion "%s" not found in elemental compositions file.\n', data.ions{iion});
    end
    
    corr_matrices{iion} = calc_corr_matrix(abundances, elements, ion_comps, ion_ind, ...
        max(cellfun(@(x) size(x, 1), data.raw(iion, :))));
end

for iion = 1:length(data.ions)
    k = strfind(data.emus{iion}, ATOM_SEP);
    natom = length(data.emus{iion}) - k(end);
    
    nrep = length(data.raw(iion, :));
    M = zeros(nrep, max(cellfun(@(x) size(x, 1), data.raw(iion, :))), length(data.t));
    if correct
        for irep = 1:nrep
            M(irep, 1:size(data.raw{iion, irep}, 1), :) ...
                = corr_matrices{iion}(1:size(data.raw{iion, irep}, 1), 1:size(data.raw{iion, irep}, 1))' ...
                \ data.raw{iion, irep};
        end
    else
        for irep = 1:nrep
            M(irep, 1:size(data.raw{iion, irep}, 1), :) = data.raw{iion, irep};
        end
    end
    
    M = M ./ repmat(sum(M, 2), [1 size(M, 2) 1]);
    M = M(:, 1:(natom + 1), :);
    
    data.mids{iion} = squeeze(mean(M, 1));
    if nrep > 1
        data.mids_std{iion} = max(squeeze(std(M, 0, 1) / sqrt(size(M, 1))), MIN_STD);   % Standard error
%         data.mids_std{iion} = max(squeeze(std(M, 0, 1)), MIN_STD);                      % Standard deviation
    else
        data.mids_std{iion} = ones(size(data.mids{iion}));
    end
        
end


fid = fopen([analysis_dir '/' SUBSTRATE_FNAME]);
C = textscan(fid, '%s%f', 'Delimiter', ',', 'Headerlines', 1);
fclose(fid);

data.src_isos = C{1};
data.src_fracs = nan(size(C{1}));
data.src_fracs(1:length(C{2})) = C{2};

end