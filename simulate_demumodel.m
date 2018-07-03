function dmids = simulate_demumodel(emumodel, v, conc0, nz_isos, ...
    nz_fracs, src_metabs, tfinal, mids, igrad)
% SIMULATE_DEMUMODEL Forward simulation of partial derivatives
%    SIMULATE_DEMUMODEL(EMUMODEL, V, CONC0, NZ_ISOS, NZ_FRACS, SRC_METABS,
%    TFINAL, IGRAD) performs forward simulation of partial derivatives of the output of EMU model EMUMODEL with
%    flux vector V, concentration vector CONC, and potential non-zero
%    isotopomers NZ_ISOS with initial fractions NZ_FRACS from time 0 to
%    TFINAL.  This function returns the partial derivative with respect to
%    component IGRAD.
%
%    The vectors V, CONC, and SRC_FRACS are column vectors, while SRC_ISOS
%    is a column cell vector.

GRAD_REL_TOL = 1e-2;
GRAD_ABS_TOL = 1e-4;

dmids = cell(length(emumodel.emus), 1);

v = reshape(v, [], 1);
conc0 = reshape(conc0, [], 1);

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

nrxn = length(v);
nmetab = length(conc0);
nnz_frac = length(nz_fracs);

for iblock = nblock:-1:1
    emu_inds = find(emu_block == iblock);
    emus_block = emus(emu_inds);
    nemus_block = length(emus_block);
    
    emu_components = strsplit(emus_block{1}, ' + ');
    if nemus_block == 1 && length(emu_components) > 1
        % EMU defined by convolution
        component_inds = cellfun(@(x) find(strcmp(x, emus)), emu_components);
        component_mids = mids(component_inds)';
        
        ncomp = length(component_inds);
        if igrad <= nrxn
            for icomp = 1:length(component_inds)
                dmids{emu_inds(1)}{1}(icomp, :) = [ dmids{component_inds(icomp)}{1} cellfun(@(x) x{1}, component_mids(setdiff(1:ncomp, icomp))) ];
                dmids{emu_inds(1)}{2}(icomp, :) = [ dmids{component_inds(icomp)}{2} cellfun(@(x) x{2}, component_mids(setdiff(1:ncomp, icomp))) ];
            end
        elseif igrad <= nrxn + nmetab
            for icomp = 1:length(component_inds)
                dmids{emu_inds(1)}{1}(icomp, :) = [ dmids{component_inds(icomp)}{1} cellfun(@(x) x{1}, component_mids(setdiff(1:ncomp, icomp))) ];
                dmids{emu_inds(1)}{2}(icomp, :) = [ dmids{component_inds(icomp)}{2} cellfun(@(x) x{2}, component_mids(setdiff(1:ncomp, icomp))) ];
            end
        else
            for icomp = 1:length(component_inds)
                dmids{emu_inds(1)}{1}(icomp, :) = [ dmids{component_inds(icomp)}{1} cellfun(@(x) x{1}, component_mids(setdiff(1:ncomp, icomp))) ];
                dmids{emu_inds(1)}{2}(icomp, :) = [ dmids{component_inds(icomp)}{2} cellfun(@(x) x{2}, component_mids(setdiff(1:ncomp, icomp))) ];
            end
        end        
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
        INSTmodel.C_diag = conc0(conc_inds);
        INSTmodel.Y = mids(emu_block ~= iblock);
        
        options = odeset('RelTol', GRAD_REL_TOL, 'AbsTol', GRAD_ABS_TOL, 'Jacobian', @INSTfjac_block);
        
        y0 = zeros(sum(emu_natoms), 1);
        INSTmodel.X = mids(emu_block == iblock);
        
        if igrad <= nrxn
            irxn = igrad;
            dA = zeros(sum(emu_natoms));
            dB = zeros(sum(emu_natoms), nnz(emu_block ~= iblock));
            
            for i = 1:nemus_block
                inds = sum(emu_natoms(1:(i - 1))) + 1;
                
                if ~ismember(emu_metabs(i), src_metabs)
                    b = ismember(A_emu_i, emu_inds) & A_emu_j == emu_inds(i);
                    in_emu_inds = A_emu_i(b);
                    coeffs = A_emu_s(b);
                    v_in_inds = A_emu_v(b);
                    dv_ins = max(coeffs .* (v_in_inds == irxn), 0);
                    
                    for j = 1:length(dv_ins)
                        j_inds = sum(emu_natoms(1:(find(emu_inds == in_emu_inds(j)) - 1))) + (1:emu_natoms(emu_inds == in_emu_inds(j)));
                        dA(inds, j_inds) = dA(inds, j_inds) + diag(dv_ins(j) * ones(length(inds), 1));
                    end
                    
                    % Subtract sum of fluxes consuming EMUs of block i
                    dA(inds, inds) = diag(sum(min(S(strcmp(emu_metabs(i), metabs), irxn), 0)) * ones(length(inds), 1));
                    
                    b = ~ismember(A_emu_i, emu_inds) & A_emu_j == emu_inds(i);
                    in_emu_inds = A_emu_i(b);
                    coeffs = A_emu_s(b);
                    v_in_inds = A_emu_v(b);
                    v_ins = max(coeffs .* (v_in_inds == irxn), 0);
                    
                    x = sparse(ones(size(in_emu_inds)), in_emu_inds, v_ins, 1, length(emus));
                    dB(inds(1), :) = x(emu_block ~= iblock);
                end
            end
            
            INSTmodel.dA = dA;
            INSTmodel.dB = dB;
            INSTmodel.dY = dmids(emu_block ~= iblock);
            
            sol = ode23t(@INSTdvode_block, [0 tfinal], y0, options);
            for i = 1:nemus_block
                dmids{emu_inds(i)}{1} = sol;
                dmids{emu_inds(i)}{2}{1} = sum(emu_natoms(1:(i - 1))) + (1:emu_natoms(i));
            end
        elseif igrad <= nrxn + nmetab
            imetab = igrad - nrxn;
            INSTmodel.dC_diag = double(conc_inds == imetab);
            INSTmodel.dY = dmids(emu_block ~= iblock);
            
            sol = ode23t(@INSTdxode_block, [0 tfinal], y0, options);
            for i = 1:nemus_block
                dmids{emu_inds(i)}{1} = sol;
                dmids{emu_inds(i)}{2}{1} = sum(emu_natoms(1:(i - 1))) + (1:emu_natoms(i));
            end
        else
            inz_frac = igrad - nrxn - nmetab;
            bnz_frac = (1:nnz_frac) == inz_frac;
            
            y0 = zeros(sum(emu_natoms), 1);
            for i = 1:nemus_block
                inds = sum(emu_natoms(1:(i - 1))) + (1:emu_natoms(i));
                
                if ismember(emu_metabs(i), nz_metabs)
                    fracs = bnz_frac(strcmp(emu_metabs(i), nz_metabs));
                    isos = nz_isos(strcmp(emu_metabs(i), nz_metabs));
                    isos_mat = cell2mat(isos);
                    
                    for j = 1:emu_natoms(i)
                        y0(inds(j)) = sum(fracs(sum(isos_mat(:, emu_atom_nums{i}), 2) == j));
                    end
                end
            end
            
            INSTmodel.dY = dmids(emu_block ~= iblock);
            
            sol = ode23t(@INSTdfode_block, [0 tfinal], y0, options);
            for i = 1:nemus_block
                dmids{emu_inds(i)}{1} = sol;
                dmids{emu_inds(i)}{2}{1} = sum(emu_natoms(1:(i - 1))) + (1:emu_natoms(i));
            end
        end
    end
end

    function J = INSTfjac_block(t, y)
        J = INSTfjac(t, y, INSTmodel);
    end

    function dydt = INSTdvode_block(t, y)
        dydt = INSTdvode_mex(t, y, INSTmodel);
%         dydt2 = INSTdvode(t, y, INSTmodel);
%         if sum(abs(dydt - dydt2)) > 1e-6
%             disp(dydt);
%             disp(dydt2);
%             save;
%         end
    end

    function dydt = INSTdxode_block(t, y)
        dydt = INSTdxode_mex(t, y, INSTmodel);
    end

    function dydt = INSTdfode_block(t, y)
        dydt = INSTdfode_mex(t, y, INSTmodel);
    end
end

