function [headers,units] = SplitHeaders(cellArray)
% ask siroos for the docs
    for i = 1 : length(cellArray)
        header = cellArray{:,i};
        [header,unit] = minor_split(header);
        headers{:,i}  = header;
        units{:,i}    = unit;
    end
end

function [header,unit] = minor_split(str)     
    header = lower(str);
    idx = regexp([' ' header],'(?<=\s+)\S','start') - 1;
    unit = str(idx(end):end);
    idx(end) = idx(end) - 2;
    header = strcat(str(idx(1):idx(end))); 
    header = regexprep(header,'(\<[a-z])','${upper($1)}');
end