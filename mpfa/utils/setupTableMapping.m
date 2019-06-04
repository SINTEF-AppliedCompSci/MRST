function [map, tbl, map1, map2] = setupTableMapping(tbl1, tbl2, crossfields, varargin)

% TODO: not supported situtation, the same fieldname is in "other fields" and
% appears as a replacement field name (see below)
    opt = struct('duplicate', [], 'fastunstable', false);
    opt = merge_options(opt, varargin{:});
    
    if ~isempty(opt.duplicate)
        dups = opt.duplicate;
        ndups = numel(dups);
        for i = 1 : ndups
            dup = dups{i};
            
            fd = dup{1};
            fd1 = dup{2}{1};            
            fd2 = dup{2}{2};
            
            ind1 = tbl1.(fd);
            ind2 = tbl2.(fd);
            
            tbl1 = rmfield(tbl1, fd);
            tbl2 = rmfield(tbl2, fd);

            tbl1.(fd1) = ind1;
            tbl2.(fd2) = ind2;
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
    
    fds1 = cell(nfds, 1);
    fds2 = cell(nfds, 1);
    for ifield = 1 : nfds
        fd = fds{ifield};
        if iscell(fd)
            % handle case with change in replacement field name.
            fds1{ifield} = fd{1};
            fds2{ifield} = fd{2};
        else
            fds1{ifield} = fd;            
            fds2{ifield} = fd;                        
        end
    end
    
    inds1   = cell(nfds, 1);
    inds2   = cell(nfds, 1);
    maxinds = cell(nfds, 1);
    prodmaxinds    = cell(nfds, 1);
    prodmaxinds{1} = 1;
    for ifield = 1 : nfds
        fieldname1 = fds1{ifield};
        fieldname2 = fds2{ifield};
        if opt.fastunstable
            inds1{ifield} = tbl1.(fieldname1);
            inds2{ifield} = tbl2.(fieldname2);
        else
            inds1{ifield} = uint64(tbl1.(fieldname1));
            inds2{ifield} = uint64(tbl2.(fieldname2));
        end
        maxinds{ifield} = max(max(inds1{ifield}), max(inds2{ifield})) + 1;
        if ifield > 1
            prodmaxinds{ifield} = prodmaxinds{ifield - 1}*maxinds{ifield - 1};
        end
    end
    
    n1 = tbl1.num; 
    n2 = tbl2.num;
    if opt.fastunstable
        globind1 = ones(n1, 1);
        globind2 = ones(n2, 1);
    else
        globind1 = ones(n1, 1, 'uint64');
        globind2 = ones(n2, 1, 'uint64');
    end
    
    for ifield = 1 : nfds
        globind1 = globind1 + inds1{ifield}*prodmaxinds{ifield};
        globind2 = globind2 + inds2{ifield}*prodmaxinds{ifield};
    end

    if opt.fastunstable
        n = max(max(globind1), max(globind2));
    else
        globind = [globind1; globind2];
        [c, ia, ic]= unique(globind);
        globind1 = ic(1 : n1);
        globind2 = ic(n1 + 1 : n1 + n2);
        n = numel(c);
    end
        
    imap1 = sparse(globind1, (1 : n1)', 1, n, n1);
    imap2 = sparse(globind2, (1 : n2)', 1, n, n2);    
    map = imap2'*imap1;

    if nargout > 1
        
        [ind2, ind1] = find(map);
    
        ofields1 = getOtherFields(tbl1, fds1);
        ofields2 = getOtherFields(tbl2, fds2);
    
        for ifield = 1 : nfds
            fieldname = fds1{ifield};
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
        
        tbl.num = numel(ind1);
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
