function tbl = sortTable(tbl, fds, varargin);
    opt = struct('keepAllFields', false);
    opt = merge_options(opt, varargin{:});
    
    a = convertTableToArray(tbl, fds);
    a = sortrows(a);
    
    if opt.keepAllFields
        warning('option not really tested yet!');
        ofds = fieldnames(tbl);
    
        fdsToRemove = {fds{:}, 'num'};
        for ifield = 1 : numel(fdsToRemove)
            ofds = ofds(~strcmp(ofds, fdsToRemove{ifield}));
        end
        
        fds = {fds{:}, ofds{:}};
    end
    
    tbl = convertArrayToTable(a, fds);
end
