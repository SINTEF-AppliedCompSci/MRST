function [map, tbl, map1, map2] = setupTableMapping(tbl1, tbl2, crossfields, varargin)
    
    opt = struct('duplicate', []);
    opt = merge_options(opt, varargin{:});
    
    if ~isempty(opt.duplicate)
        dups = opt.duplicate;
        ndups = numel(dups);
        for i = 1 : ndups
            dup = dups{i};
            fd = dup{1};
            fd1 = dup{2}{1};            
            fd2 = dup{2}{2};
            tbl1.(fd1) = tbl1.(fd);
            tbl1 = rmfield(tbl1, fd);
            tbl2.(fd2) = tbl2.(fd);
            tbl2 = rmfield(tbl2, fd);
        end
    end
    
    fds = crossfields; %alias
    nfds = numel(fds);
    
    if isempty(crossfields)
        % handle case when crossfield is empty
        n1 = tbl1.num;
        n2 = tbl2.num;
        map = ones(n2, n1);
    end
    
    for ifield = 1 : nfds
        fieldname = fds{ifield};
        
        ind1 = tbl1.(fieldname);
        ind2 = tbl2.(fieldname);
        
        n1 = numel(ind1);
        n2 = numel(ind2);
        n = max([ind1; ind2]);
        
        imap1 = sparse(ind1, (1 : n1)', 1, n, n1);
        imap2 = sparse(ind2, (1 : n2)', 1, n, n2);
        
        imap = imap2'*imap1;
        if ifield == 1
            map = imap;
        else
            map = map.*imap;
        end
    end

    if nargout > 1
        [ind2, ind1] = find(map);
    
        ofields1 = getOtherFields(tbl1, fds);
        ofields2 = getOtherFields(tbl2, fds);
    
        for ifield = 1 : nfds
            fieldname = fds{ifield};
            tbl.(fieldname) = tbl1.(fieldname)(ind1);
        end
        
        for ifield = 1 : numel(ofields1)
            fieldname = ofields1{ifield};
            tbl.(fieldname) = tbl1.(fieldname)(ind1);
        end

        for ifield = 1 : numel(ofields2)
            fieldname = ofields2{ifield};
            tbl.(fieldname) = tbl2.(fieldname)(ind2);
        end
        
        % get dimension of table (may not be so robust, to be check...).
        tbl.num = numel(tbl.(ofields1{1}));
    end
    
    if nargout > 2
        error('not working for moment')
        map1 = setupTableMapping(tbl1, tbl, fds);
        map2 = setupTableMapping(tbl2, tbl, fds);
    end
    
end

function ofields = getOtherFields(tbl, fields)
    ofields = fieldnames(tbl);
    if any(strcmp(ofields, 'ind'));
        ofields = ofields(~strcmp(ofields, 'ind'));
    end
    if any(strcmp(ofields, 'num'));
        ofields = ofields(~strcmp(ofields, 'num'));
    end    
    ofields = ofields(~ismember(ofields, fields));
end
