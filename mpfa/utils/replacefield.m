function newtbl = replacefield(tbl, oldfield, newfield)
    newtbl = tbl;
    newtbl.(newfield) = tbl.(oldfield);
    newtbl = rmfield(newtbl, oldfield);
end
