function cl = fullStructToCell(str)
% Convert structure to cell array consisting of interleaved field names and
% field values.  This function thus differ from the in-built Matlab function
% 'str2cell' by also including the fieldnames in the cell array.

    cl = interleave(fieldnames(str), struct2cell(str));
    
end
