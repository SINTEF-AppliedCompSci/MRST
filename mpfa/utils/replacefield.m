function newtbl = replacefield(tbl, fieldpairs)
    newtbl = tbl;
    if iscell(fieldpairs{1})
        for i = 1 : numel(fieldpairs)
            fielpair = fieldpairs{i};
            newfield = fieldpair{1};
            oldfield = fieldpair{2};
            newtbl.(newfield) = tbl.(oldfield);
            newtbl = rmfield(newtbl, oldfield);
        end
    else
        fieldpair = fieldpairs{1};
        newfield = fieldpair{1};
        oldfield = fieldpair{2};
        newtbl.(newfield) = tbl.(oldfield);
        newtbl = rmfield(newtbl, oldfield);
    end
end
