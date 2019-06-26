function tbl = addLocInd(tbl, locindname)
    tbl.(locindname) = (1 : tbl.num)';
end
