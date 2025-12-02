function [data,table] = ImportTable(file_path)
%
% DESCRIPTION: reads tabular data from kr, pc, schedule, ... tables in .txt
%              format
%
% SYNOPSIS:
%   [data,table] = ImportTable(fullFile)
%
% PARAMETERS:
%   fullFile - string path + file name (full file path) to the .txt file
%
% RETURNS:
%   data - array of the read data
%   table - table of the read data
%
% EXAMPLE:
%   [data,table] = ImportTable(".\kr_table.txt")
%
% ----------------------------------
% (c) 2020-2022
% Siroos Azizmohammadi
% Omidreza Amrollahinasab
% Montanuniversität Leoben, Austria
% Chair of Reservoir Engineering
% https://dpe.ac.at/
% ----------------------------------
%
%%
file = fopen(file_path,'r');
str = fileread(file_path); % read entire file into string
lines = strtrim(regexp( str,'(\r|\n)+','split')); % split by each line
data = [];
cols = regexp(string(lines(2)),'''(.[^'']*)''','tokens');
headers = [];
for i = 1 : length(cols)
    headers = [headers, cols{i}];
    newUnits{i} = Unit(cols{i});
end
cols = regexp(string(lines(4)),'''(.[^'']*)''','tokens');
units = [];
for i = 1 : length(cols)
    units = [units, cols{i}];
end
tableHeaders = [];
for i = 1 : length(cols)
    tableHeaders = [tableHeaders, strcat(headers(i),' [',newUnits(i),']')];
end
headerLines = 4;
for k = headerLines + 1 : length(lines)            
    line = strtrim(regexp(lines{k},'\s+','split')); % split by spaces
    while (~isempty(find(~cellfun(@isempty,line),1)) && ...
           ~startsWith(line{1},'#'))
        row = [];
        for j = 1 : length(cols)
            value = str2double(line(j));
            unit = units(j);
            row = [row, value * Convert(unit)];
        end
        data = [data;row];
        break
    end
end
table = array2table(data,'VariableNames',tableHeaders);
fclose(file);