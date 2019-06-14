function [restbl, map1, map2, reducemap] = setupTableMultOperators(tbl1, tbl2, ...
                                                                         fds1, ...
                                                                         fds2, ...
                                                                         fdscross)
    afds1 = fieldnames(tbl1);
    afds2 = fieldnames(tbl2);
    
    %% sanity checks
    % first check: we should have fds1 (fds2) equal to ofds1 (ofds2).
    ofds1 = afds1;
    fdsToRemove = {fdscross{:}, 'ind', 'num'};
    for ifield = 1 : numel(fdsToRemove)
        ofds1 = ofds1(~strcmp(ofds1, fdsToRemove{ifield}));
    end
    
    ofds2 = afds2;
    fdsToRemove = {fdscross{:}, 'ind', 'num'};
    for ifield = 1 : numel(fdsToRemove)
        ofds2 = ofds2(~strcmp(ofds2, fdsToRemove{ifield}));
    end
    
    isdiff = any(setdiff(fds1, ofds1)) | any(setdiff(fds2, ofds2));
    assert(~isdiff, 'mismatch in table matrix multipication');
    % second check: we do not support for the moment when fds1 and fds2 have
    % common field names (in this case they should be given different names,
    % which we do not do automatically).
    
    assert(nnz(intersect(ofds1, ofds2)) == 0, ['repeated names in fields are not ' ...
                        'supported']);
    
    % end of sanity checks
    [~, prodmattbl] = setupTableMapping(tbl1, tbl2, fdscross);
    map1 = setupTableMapping(tbl1, prodmattbl, fdscross);
    map2 = setupTableMapping(tbl2, prodmattbl, fdscross);
    resfds = {fds1{:}, fds2{:}};
    restbl = projTable(prodmattbl, resfds);
    reducemap = setupTableMapping(prodmattbl, restbl, resfds);
    
end
