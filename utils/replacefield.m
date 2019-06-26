function newtbl = replacefield(tbl, fieldpairs)
    newtbl = tbl;
    if iscell(fieldpairs{1})
        for i = 1 : numel(fieldpairs)
            newtbl = replacethisfield(newtbl, fieldpairs{i});
        end
    else
        fieldpair = fieldpairs;
        newtbl = replacethisfield(newtbl, fieldpair);
    end
end

function tbl = replacethisfield(tbl, fieldpair)
    oldfield = fieldpair{1};
    newfield = fieldpair{2};
    tbl.(newfield) = tbl.(oldfield);
    tbl = rmfield(tbl, oldfield);
end