function [data,table] = ImportTable(fullFile)
    file = fopen(fullFile,'r');
    str = fileread(fullFile); % read entire file into string
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
end