function emumodel = generate_emumodel(in_fname, measured_emus, bnd_metabs)
% GENERATE_EMUMODEL Generates a structure that contains the information
% about the EMU reactions taking part in the network.
%    EMUMODEL = GENERATE_EMUMODEL(IN_FNAME, MEASURED_EMUS, BND_METABS)
%    generates an EMU model from the model input file IN_FNAME for the set
%    of measured EMUs in the cell array MEASURED_EMUS and the set of 
%    boundary metabolites in the cell array BND_METABS.
%
%    The EMU model is generated using the reaction network decoupling
%    method proposed by Young et al. (Biotech. and Bioeng., 99(3):686-99,
%    2007).
%
%    Example:
%       emumodel = generate_emumodel('test.csv', {'C_1234'}, {'A', 'J'})
% Metabolite names cannot contain ATOM_SEP or COPY_SEP.

ATOM_SEP = '_';
COPY_SEP = ':';
COMPART_SEP = '.';

if nargin < 3
    bnd_metabs = {};
end

fid = fopen(in_fname);
C = textscan(fid, '%s%f', 'Delimiter', ',', 'TreatAsEmpty', {'None'});
fclose(fid);

nrxn = length(C{1});
reverse = false(nrxn, 1);
metabs = {};
imetab = 1;
atoms = {};
iatom = 1;
metab_equiv = {};

Si = []; Sj = []; Ss = [];  
Ai = []; Aj = []; As = [];  Av = [];

for irxn = 1:nrxn
    rxnstring = C{1}{irxn};
    
    try
        if ~isempty(strfind(rxnstring, '<==>'))
            eqsym = '<==>';
            reverse(irxn) = true;
        elseif ~isempty(strfind(rxnstring, '-->'))
            eqsym = '-->';
            reverse(irxn) = false;
        elseif ~isempty(strfind(rxnstring, '<=>'))
            eqsym = '<=>';
            reverse(irxn) = true;
        elseif ~isempty(strfind(rxnstring, '<->'))
            eqsym = '<->';
            reverse(irxn) = true;
        elseif ~isempty(strfind(rxnstring, '=>'))
            eqsym = '=>';
            reverse(irxn) = false;
        elseif ~isempty(strfind(rxnstring, '->'))
            eqsym = '->';
            reverse(irxn) = false;
        elseif ~isempty(strfind(rxnstring, '<-'))
            eqsym = '<-';
            reverse(irxn) = false;
        elseif ~isempty(strfind(rxnstring, '='))
            eqsym = '=';
            reverse(irxn) = true;
        end
        
        labeled_stoich = Inf;
        k = strfind(rxnstring, eqsym);
        lhs = strtrim(rxnstring(1:(k(1) - 1)));
        lhs_trans = '';
        lhs_atoms = {};
        if ~isempty(lhs)
            terms = strsplit(lhs, ' + ');
            terms = strtrim(terms);
            n = regexp(terms, '(?<stoich>\(?[0-9.]+\)?\s+|)(?<metab>\S+)(?<trans>\s+\([^\)]+\)|)', 'names');
            for iterm = 1:length(n)
                if ~isempty(n{iterm}.stoich)
                    s = regexp(n{iterm}.stoich, '\(?([0-9.]+)\)?', 'tokens');
                    stoich = str2double(s{1});
                else
                    stoich = 1;
                end
                metab = n{iterm}.metab;
                
                i = find(strcmp(metab, metabs));
                if isempty(i)
                    if any(metab == ATOM_SEP) || any(metab == COPY_SEP)
                        error('Invalid metabolite "%s". Metabolites cannot contain "%c" or "%c".', ...
                            metab, ATOM_SEP, COPY_SEP);
                    end
                    
                    i = imetab;
                    metabs{imetab} = metab;
                    metab_equivs{imetab} = [];
                    imetab = imetab + 1;
                end
                Si = [Si; i];
                Sj = [Sj; irxn];
                Ss = [Ss; -stoich];
                
                if ~isempty(n{iterm}.trans)
                    labeled_stoich = min(labeled_stoich, stoich);
                    
                    s = regexp(n{iterm}.trans, '\(([^\)]+)\)', 'tokens');
                    trans_terms = strsplit(s{1}{1}, ' = ');
                    lhs_trans = strcat(lhs_trans, trans_terms{1});
                    natoms = length(trans_terms{1});
                    lhs_atoms = [lhs_atoms; cellstr([repmat(strcat(metab, ATOM_SEP), natoms, 1) char('0' + (1:natoms)')])];
                    
                    for j = 2:length(trans_terms)
                        metab_equivs{i}{j - 1} = arrayfun(@(x) find(x == trans_terms{j}), trans_terms{1});
                    end
                end
            end
        end
        
        rhs = strtrim(rxnstring((k(1) + length(eqsym)):end));
        rhs_trans = '';
        rhs_atoms = {};
        if ~isempty(rhs)
            terms = strsplit(rhs, ' + ');
            terms = strtrim(terms);
            n = regexp(terms, '(?<stoich>\(?[0-9.]+\)?\s+|)(?<metab>\S+)(?<trans>\s+\([^\)]+\)|)', 'names');
            for iterm = 1:length(n)
                if ~isempty(n{iterm}.stoich)
                    s = regexp(n{iterm}.stoich, '\(?([0-9.]+)\)?', 'tokens');
                    stoich = str2double(s{1});
                else
                    stoich = 1;
                end
                metab = n{iterm}.metab;
                
                i = find(strcmp(metab, metabs));
                if isempty(i)
                    if any(metab == ATOM_SEP) || any(metab == COPY_SEP)
                        error('Invalid metabolite "%s". Metabolites cannot contain "%c" or "%c".', ...
                            metab, ATOM_SEP, COPY_SEP);
                    end
                    
                    i = imetab;
                    metabs{imetab} = metab;
                    metab_equivs{imetab} = [];
                    imetab = imetab + 1;
                end
                Si = [Si; i];
                Sj = [Sj; irxn];
                Ss = [Ss; stoich];
                
                if ~isempty(n{iterm}.trans)
                    labeled_stoich = min(labeled_stoich, stoich);
                    
                    s = regexp(n{iterm}.trans, '\(([^\)]+)\)', 'tokens');
                    trans_terms = strsplit(s{1}{1}, ' = ');
                    rhs_trans = strcat(rhs_trans, trans_terms{1});
                    natoms = length(trans_terms{1});
                    rhs_atoms = [rhs_atoms; cellstr([repmat(strcat(metab, ATOM_SEP), natoms, 1) char('0' + (1:natoms)')])];
                    
                    for j = 2:length(trans_terms)
                        metab_equivs{i}{j - 1} = arrayfun(@(x) find(x == trans_terms{j}), trans_terms{1});
                    end
                end
            end
        end
        
        [c, ia, ic] = unique(lhs_atoms);
        if length(ia) ~= length(ic)
            for atom = reshape(c, 1, [])
                batom = strcmp(atom, lhs_atoms);
                lhs_atoms(batom) = strcat(lhs_atoms(batom), cellstr(char(COPY_SEP * ones(nnz(batom), 1))), cellstr(char('0' + (1:nnz(batom)))'));
            end
        end
        
        [c, ia, ic] = unique(rhs_atoms);
        if length(ia) ~= length(ic)
            for atom = reshape(c, 1, [])
                batom = strcmp(atom, rhs_atoms);
                rhs_atoms(batom) = strcat(rhs_atoms(batom), cellstr(char(COPY_SEP * ones(nnz(batom), 1))), cellstr(char('0' + (1:nnz(batom)))'));
            end
        end
        
        lr_atoms = [lhs_atoms; rhs_atoms];
        for jatom = 1:length(lr_atoms)
            i = find(strcmp(lr_atoms{jatom}, atoms));
            if isempty(i)
                i = iatom;
                atoms{iatom} = lr_atoms{jatom};
                iatom = iatom + 1;
            end
        end
        
        ilhs_atoms = cellfun(@(x) find(strcmp(x, atoms)), lhs_atoms);
        irhs_atoms = cellfun(@(x) find(strcmp(x, atoms)), rhs_atoms);
        
        if ~isempty(lhs_atoms) && ~isempty(rhs_atoms)
            Ai = [Ai; ilhs_atoms];
            Aj = [Aj; irhs_atoms(arrayfun(@(x) find(x == rhs_trans), lhs_trans))];
            Av = [Av; irxn * ones(length(ilhs_atoms), 1)];
            As = [As; labeled_stoich * ones(length(ilhs_atoms), 1)];
            
            if reverse(irxn)
                Ai = [Ai; irhs_atoms(arrayfun(@(x) find(x == rhs_trans), lhs_trans))];
                Aj = [Aj; ilhs_atoms];
                Av = [Av; -irxn * ones(length(ilhs_atoms), 1)];
                As = [As; labeled_stoich * ones(length(ilhs_atoms), 1)];
            end
        end
    catch
        error('Error parsing reaction "%s".', rxnstring);
    end
end

S = sparse(Si, Sj, Ss, length(metabs), nrxn);               
vfixed = C{2};
nfixed = nnz(~isnan(vfixed));
Afixed = sparse(1:nfixed, find(~isnan(vfixed)), ones(nfixed, 1), nfixed, nrxn);
bfixed = vfixed(~isnan(vfixed));

% Convert reversible reactions to two irreversible reactions
S = [S -S(:, reverse)];
Afixed = [Afixed -Afixed(:, reverse)];
irev_rxn = 1;
for irxn = find(reverse)'
    Av(Av == -irxn) = nrxn + irev_rxn;
    irev_rxn = irev_rxn + 1;
end

% Expand measured_emus to compartment instances
tmp_emus = {};
for iemu = 1:length(measured_emus)
    emu = measured_emus{iemu};
    k = strfind(emu, ATOM_SEP);
    if isempty(k)
        error('Invalid measured EMU "%s".  Must contain "%c".', ...
            emu, ATOM_SEP);
    end
    
    emu_metab = emu(1:(k(end) - 1));
    atom_nums = emu((k(end) + 1):end);
    for metab = metabs(strncmp(emu_metab, metabs, length(emu_metab)))
        if length(metab{1}) == length(emu_metab) || ...
            metab{1}(length(emu_metab) + 1) == COMPART_SEP
            tmp_emus(end + 1) = strcat(metab, ATOM_SEP, atom_nums);
        end
    end
end
measured_emus = tmp_emus;

% Generate EMUs
emu_stack = {};
for iemu = 1:length(measured_emus)
    emu = measured_emus{iemu};
    k = strfind(emu, ATOM_SEP);
    emu_size = length(emu) - k(end);
    if length(emu_stack) < emu_size || isempty(emu_stack(emu_size))
        emu_stack{emu_size} = { emu };
    else
        emu_stack{emu_size}{end + 1} = emu;
    end
end

emus = measured_emus;
iemu = length(emus) + 1;
A_emu_i = []; A_emu_j = []; A_emu_s = []; A_emu_v = []; 
emus_visited = {};
while sum(cellfun(@length, emu_stack)) > 0  
    emu_size = find(cellfun(@length, emu_stack) > 0, 1, 'last');
    emu = emu_stack{emu_size}{1};
    
    emu_components = strsplit(emu, ' + ');    
    if length(emu_components) > 1
        % EMU is a convolution
        for icomponent = 1:length(emu_components)
            in_emu = emu_components{icomponent};
            
            i = find(strcmp(in_emu, emus));
            if isempty(i)
                i = iemu;
                emus{iemu} = in_emu;
                iemu = iemu + 1;

                k = strfind(in_emu, ATOM_SEP);
                in_metab = in_emu(1:(k(end) - 1));
                if ~ismember(in_metab, bnd_metabs) 
                    % Add in_emu to emu_stack
                    in_emu_size = length(in_emu) - k(end);
                    if length(emu_stack) < in_emu_size || isempty(emu_stack(in_emu_size))
                        emu_stack{in_emu_size} = { in_emu };
                    else
                        emu_stack{in_emu_size}{end + 1} = in_emu;
                    end
                end
            end
            
            A_emu_i = [A_emu_i; i];
            A_emu_j = [A_emu_j; find(strcmp(emu, emus))];
            A_emu_s = [A_emu_s; 1];
            A_emu_v = [A_emu_v; 0];
        end
    else
        k = strfind(emu, ATOM_SEP);
        metab = emu(1:(k(end) - 1));  
        metab_equiv = metab_equivs{strcmp(metab, metabs)};
        for iequiv = 1:(length(metab_equiv) + 1)
            atom_set = emu((k(end) + 1):end);
            if iequiv > 1 
                batoms = false(length(metab_equiv{iequiv - 1}), 1);
                batoms(atom_set - '0') = true;

                atom_set = char(find(batoms(metab_equiv{iequiv - 1}))' + '0');
            end
            emu_atoms = cellstr([ repmat([ metab ATOM_SEP ], length(atom_set), 1) atom_set' ]);
            emu_atoms_inds = cellfun(@(x) find(strcmp(x, atoms)), emu_atoms);

            while ~isempty(emu_atoms_inds)
                irxns = unique(Av(ismember(Aj, emu_atoms_inds)))';  % indices of incoming reactions to EMU
                for irxn = irxns
                    i = Ai(ismember(Aj, emu_atoms_inds) & Av == irxn);
                    in_atoms = atoms(i);
                    stoich = As(ismember(Aj, emu_atoms_inds) & Av == irxn);
                    stoich = stoich(1);

                    n = regexp(in_atoms, '(?<metab>\S+)_(?<atom_num>\S)(?<rep_num>.\S)?', 'names');
                    in_metabs = cellfun(@(x) x.metab, n, 'UniformOutput', false);
                    in_atom_nums = cellfun(@(x) x.atom_num, n);
                    in_rep_nums = cellfun(@(x) x.rep_num, n, 'UniformOutput', false);
                    in_metabs_rep = strcat(in_metabs, in_rep_nums);

                    unique_in_metabs_rep = unique(in_metabs_rep);
                    in_emus = cell(length(unique_in_metabs_rep), 1);
                    for imetab = 1:length(unique_in_metabs_rep)
                        in_metab_rep = unique_in_metabs_rep{imetab};                
                        in_atom_set = sort(in_atom_nums(strcmp(in_metab_rep, in_metabs_rep)));
                        
                        if length(in_metab_rep) <= 1 || in_metab_rep(end - 1) ~= COPY_SEP
                            in_metab = in_metab_rep;
                        else
                            in_metab = in_metab_rep(1:(end - 2));
                        end
                        in_metab_equiv = metab_equivs{strcmp(in_metab, metabs)};
                        for iin_equiv = 1:(length(in_metab_equiv) + 1)
                            if iin_equiv > 1
                                batoms = false(length(in_metab_equiv{iin_equiv - 1}), 1);
                                batoms(in_atom_set - '0') = true;

                                in_atom_set = char(find(batoms(in_metab_equiv{iin_equiv - 1}))' + '0');
                            end
                            
                            in_emu = [ in_metab ATOM_SEP in_atom_set ];
                            
                            if any(strcmp(in_emu, emus))
                                break;
                            end
                        end

                        in_emus{imetab} = in_emu;               
                    end

                    full_in_emu = in_emus{1};
                    for i = 2:length(in_emus)
                        full_in_emu = [ full_in_emu ' + ' in_emus{i} ];
                    end

                    i = find(strcmp(full_in_emu, emus));
                    if isempty(i)
                        i = iemu;
                        emus{iemu} = full_in_emu;
                        iemu = iemu + 1;

                        if ~all(ismember(unique_in_metabs_rep, bnd_metabs))
                            % Add in_emu to emu_stack
                            in_emu_size = nnz(strcmp(in_metab, in_metabs));
                            if length(emu_stack) < in_emu_size || isempty(emu_stack(in_emu_size))
                                emu_stack{in_emu_size} = { full_in_emu };
                            else
                                emu_stack{in_emu_size}{end + 1} = full_in_emu;
                            end
                        end
                    end

                    A_emu_i = [A_emu_i; i];
                    A_emu_j = [A_emu_j; find(strcmp(emu, emus))];
                    A_emu_s = [A_emu_s; stoich / (length(metab_equiv) + 1)];
                    A_emu_v = [A_emu_v; irxn];     
                end
                
                if emu_atoms{1}((end - 1)) ~= COPY_SEP
                    emu_atoms = strcat(emu_atoms, cellstr(repmat([ COPY_SEP '1' ], length(emu_atoms), 1)));
                else
                    emu_atoms = cellfun(@(x) [x(1:(end - 1)) x(end) + 1], emu_atoms, 'UniformOutput', false);
                end
                emu_atoms_inds = cell2mat(cellfun(@(x) find(strcmp(x, atoms)), emu_atoms, 'UniformOutput', false));
            end
        end
    end
    
    emu_stack{emu_size} = setdiff(emu_stack{emu_size}, emu);
end

A_emu = sparse(A_emu_i, A_emu_j, A_emu_s, length(emus), length(emus));
[nblock, emu_block] = graphconncomp(A_emu);

emumodel.emus = reshape(emus, [], 1);           
emumodel.A_emu_i = A_emu_i;     
emumodel.A_emu_j = A_emu_j;
emumodel.A_emu_s = A_emu_s;
emumodel.A_emu_v = A_emu_v;     
emumodel.nblock = nblock;       
emumodel.emu_block = emu_block;

emumodel.S = S;                 
emumodel.metabs = metabs;      
emumodel.reverse = reverse;

emumodel.Afixed = Afixed;
emumodel.bfixed = bfixed;

end

