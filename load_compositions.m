function ion_comps = load_compositions(filename)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

ion_comps = struct();

try
    data = importdata(filename);
catch
    error('Unable to load file "%s"\n', filename);
end

header_rows = 1; %(There is a header row...)

nions = size(data,1)-header_rows;
ion_comps.ions = cell(nions,1);
ion_comps.elements = cell(nions,1);
ion_comps.numbers = cell(nions,1);

for i = 1:(size(data,1)-header_rows)
    C = strsplit(data{i+header_rows,1},',');
    
    [ele_list,~,ele_end] = regexp(C{2},['[','A':'Z','][','a':'z',']?'],'match');
    
    [num,num_start] = regexp(C{2},'\d+','match');
    num_list = ones(1,length(ele_list));
    index = ismember(ele_end+1, num_start);
    num_list(index) = cellfun(@str2num, num);
    
    ion_comps.ions{i} = C{1};
    ion_comps.elements{i} = ele_list;
    ion_comps.numbers{i} = num_list;
end

end

