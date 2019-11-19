function y = convertToColumn(y)
    if size(y,1) < size(y,2)
        y = y';
    end
end