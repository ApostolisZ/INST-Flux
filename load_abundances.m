function [elements, abundances] = load_abundances(filename)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

try
    data = importdata(filename);
catch
    error('Unable to load file "%s"\n', filename);
end

elements = data.textdata;
abundances = cell(size(elements)); % textdata and rowheaders are the same
                                   % and stores the first textual column
                                   % from the file
                                   
for i = 1:size(elements,1)
    abundances{i} = nonzeros(data.data(i,2:end));
    %disp(abundances{i})
end

end

