function [ corr_matrix ] = calc_corr_matrix( abundances, elements, composition, iion, matrix_size )
% [CORR_MATRIX] = CALC_CORR_MATRIX( ABUNDANCES, ELEMENTS, COMPOSITIONS, IION, MATRIX_SIZE )
% Calculates correction matrices, to adjust labeling for natural
% abundance
%   CORR_MATRIX is a matrix containing the natural labeling abundance of a
%   metabolite with (MATRIX_SIZE-1) atoms.
%   ABUNDANCES is the natural labeling abundances for the chemical elements
%   present in the metabolite
%   ELEMENTS contains the chemical elements present in the network
%   COMPOSITION contains the chemical elements present in the metabolite
%   IION is the index of the labeled EMU
%   MATRIX_SIZE is the number of labeling states of the labeling states

classical = 0;

corr_matrix = eye(matrix_size);

if classical
    for element = 1:size(composition.elements{iion},1)
        index = strcmp(elements,composition.elements{iion}(element));
        if any(index)
            element_corr_matrix = zeros(matrix_size);
            matrix_row = corr_matrix_row(abundances{index}',composition.numbers{iion}(element),matrix_size);
            for i = 1:matrix_size %rows
                last_pos = matrix_size + 1 - i;
                element_corr_matrix(i,i:end) = matrix_row(1:last_pos);
            end
            corr_matrix = corr_matrix * element_corr_matrix;
        end
    end
else
    for element = 1:size(composition.elements{iion},2)
        index = strcmp(elements,composition.elements{iion}(element));
        if any(index)
            element_corr_matrix = zeros(matrix_size);
            if strcmp(composition.elements{iion}(element),'C')
                for i = 1:matrix_size %rows
                    matrix_row = corr_matrix_row(abundances{index}',composition.numbers{iion}(element)-i+1,matrix_size);
                    last_pos = matrix_size + 1 - i;
                    element_corr_matrix(i,i:end) = matrix_row(1:last_pos);
                end
            else
                matrix_row = corr_matrix_row(abundances{index}',composition.numbers{iion}(element),matrix_size);
                for i = 1:matrix_size %rows
                    last_pos = matrix_size + 1 - i;
                    element_corr_matrix(i,i:end) = matrix_row(1:last_pos);
                end
            end
            corr_matrix = corr_matrix * element_corr_matrix;
        end
    end
end

end

