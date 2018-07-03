function mids = simulate_emumodel(emumodel, v, conc, nz_isos, nz_fracs, src_metabs, tfinal)
% SIMULATE_EMUMODEL Forward simulation of EMU model
%    SIMULATE_EMUMODEL(EMUMODEL, V, CONC, NZ_ISOS, NZ_FRACS, SRC_METABS,
%    TFINAL) performs forward simulation of EMU model EMUMODEL with flux
%    vector V, concentration vector CONC, and potential non-zero
%    isotopomers NZ_ISOS with initial fractions NZ_FRACS from time 0 to
%    TFINAL.
%
%    The vectors V, CONC, and SRC_FRACS are column vectors, while SRC_ISOS
%    is a column cell vector.

REL_TOL = 1e-3;
ABS_TOL = 1e-5;

mids = cell(length(emumodel.emus), 1);

v = reshape(v, [], 1);
conc = reshape(conc, [], 1);

emus = emumodel.emus;
A_emu_i = emumodel.A_emu_i;
A_emu_j = emumodel.A_emu_j;
A_emu_s = emumodel.A_emu_s;
A_emu_v = emumodel.A_emu_v;
nblock = emumodel.nblock;
emu_block = emumodel.emu_block;

S = emumodel.S;
metabs = emumodel.metabs;

n = regexp(nz_isos, '(?<metab>\S+)_(?<iso>\S+)', 'names');
nz_metabs = cellfun(@(x) x.metab, n, 'UniformOutput', false);
nz_isos = cellfun(@(x) x.iso - '0', n, 'UniformOutput', false);

for iblock = nblock:-1:1
    emu_inds = find(emu_block == iblock);
    emus_block = emus(emu_inds);
    nemus_block = length(emus_block);
    
    emu_components = strsplit(emus_block{1}, ' + ');
    if nemus_block == 1 && length(emu_components) > 1
        % EMU defined by convolution
        component_inds = cellfun(@(x) find(strcmp(x, emus)), emu_components);
        component_mids = mids(component_inds)';
        
        mids{emu_inds(1)}{1} = cellfun(@(x) x{1}, component_mids);
        mids{emu_inds(1)}{2} = cellfun(@(x) x{2}, component_mids);
    else
        n = regexp(emus_block, '(?<metab>\S+)_(?<atom_nums>\S+)', 'names');
        emu_metabs = cellfun(@(x) x.metab, n, 'UniformOutput', false);    
        emu_atom_nums = cellfun(@(x) x.atom_nums - '0', n, 'UniformOutput', false);
        emu_natoms = cellfun(@(x) length(x.atom_nums), n);
        INSTmodel.emu_natoms = emu_natoms;

        conc_inds = zeros(sum(emu_natoms), 1);
        y0 = zeros(sum(emu_natoms), 1);
        A = zeros(sum(emu_natoms));
        B = zeros(sum(emu_natoms), nnz(emu_block ~= iblock));
        for i = 1:nemus_block
            inds = sum(emu_natoms(1:(i - 1))) + (1:emu_natoms(i));
            conc_inds(inds) = find(strcmp(emu_metabs(i), metabs));

            if ismember(emu_metabs(i), nz_metabs)
                fracs = nz_fracs(strcmp(emu_metabs(i), nz_metabs));
                isos = nz_isos(strcmp(emu_metabs(i), nz_metabs));
                isos_mat = cell2mat(isos);

                for j = 1:emu_natoms(i)
                    y0(inds(j)) = sum(fracs(sum(isos_mat(:, emu_atom_nums{i}), 2) == j));
                end
            end
            
            if ~ismember(emu_metabs(i), src_metabs)    
                b = ismember(A_emu_i, emu_inds) & A_emu_j == emu_inds(i);
                in_emu_inds = A_emu_i(b);
                coeffs = A_emu_s(b);
                v_in_inds = A_emu_v(b);
                v_ins = max(coeffs .* v(v_in_inds), 0);

                for j = 1:length(v_ins)
                    j_inds = sum(emu_natoms(1:(find(emu_inds == in_emu_inds(j)) - 1))) + (1:emu_natoms(emu_inds == in_emu_inds(j)));
                    A(inds, j_inds) = A(inds, j_inds) + diag(v_ins(j) * ones(length(inds), 1));
                end

                % Subtract sum of fluxes consuming EMUs of block i
                A(inds, inds) = diag(sum(min(S(strcmp(emu_metabs(i), metabs), :) .* v', 0)) * ones(length(inds), 1));
                
                b = ~ismember(A_emu_i, emu_inds) & A_emu_j == emu_inds(i);
                in_emu_inds = A_emu_i(b);
                coeffs = A_emu_s(b);
                v_in_inds = A_emu_v(b);
                v_ins = max(coeffs .* v(v_in_inds), 0);
                
                x = sparse(ones(size(in_emu_inds)), in_emu_inds, v_ins, 1, length(emus));
                B(inds(1), :) = x(emu_block ~= iblock);
            end
        end
        
        INSTmodel.A = A;
        INSTmodel.B = B;
        INSTmodel.C_diag = conc(conc_inds);
        INSTmodel.Y = mids(emu_block ~= iblock);
        
        options = odeset('RelTol', REL_TOL, 'AbsTol', ABS_TOL, 'Jacobian', @INSTfjac_block);
        sol = ode23t(@INSTode_block, [0 tfinal], y0, options);
        
        %         disp(emus_block);
        %         plot(sol.x, sol.y, 'o-');
        %         pause;
        for i = 1:nemus_block
            mids{emu_inds(i)}{1} = sol;
            mids{emu_inds(i)}{2}{1} = sum(emu_natoms(1:(i - 1))) + (1:emu_natoms(i));
        end
    end
end

    function dydt = INSTode_block(t, y)
        dydt = INSTode_mex(t, y, INSTmodel);
    end

    function J = INSTfjac_block(t, y)
        J = INSTfjac(t, y, INSTmodel);
    end
end

