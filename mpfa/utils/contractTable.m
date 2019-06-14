function [prodAB, restbl] = multTable(Acell, Bcell, fds)
    
    A    = Acell{1};
    tblA = Acell{2};
    B    = Bcell{1};
    tblB = Bcell{2};
    
    fds1 = fds{1};
    fds2 = fds{2};
    fdscross = fds{3};
    
    afds1 = fieldnames(tblA);
    afds2 = fieldnames(tblB);
    
    %% sanity checks
    % first check: we should have fds1 (fds2) included ofds1 (ofds2).
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
    
    isdiff = (numel(setdiff(fds1, ofds1)) > 0) | (numel(setdiff(fds2, ofds2)) > 0);
    assert(~isdiff, 'mismatch in table matrix multipication');
    % second check: we do not support for the moment when fds1 and fds2 have
    % common field names (in this case they should be given different names,
    % which we do not do automatically).
    
    assert(numel(intersect(ofds1, ofds2)) == 0, ['repeated names in fields are not ' ...
                        'supported']);
    
    % end of sanity checks
    
    [~, prodmattbl] = setupTableMapping(tblA, tblB, fdscross);
    map1 = setupTableMapping(tblA, prodmattbl, {ofds1{:}, fdscross{:}});
    map2 = setupTableMapping(tblB, prodmattbl, {ofds2{:}, fdscross{:}});
    resfds = {fds1{:}, fds2{:}};
    restbl = projTable(prodmattbl, resfds);
    reducemap = setupTableMapping(prodmattbl, restbl, resfds);
    
    prodAB = reducemap*((map1*A).*(map2*B));
end
