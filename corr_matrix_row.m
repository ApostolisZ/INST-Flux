function [ matrix_row ] = corr_matrix_row( abundances, number, matrix_size )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
% Image Processing Toolbox is required to implement padarray dealing with input
% type of 'double'

matrix_row = 1;
if number > 0
    for i = 1:number
        matrix_row = conv(matrix_row,abundances);
    end
end

if matrix_size > length(matrix_row)
    matrix_row = padarray(matrix_row,[0 (matrix_size - length(matrix_row))], 'post');
else
    matrix_row = matrix_row(1:matrix_size);
end

end

