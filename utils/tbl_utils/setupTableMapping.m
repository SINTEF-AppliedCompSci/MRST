function [map, tbl] = setupTableMapping(tbl1, tbl2, crossfields, varargin)
%
%
% SYNOPSIS:
%   function [map, tbl, map1, map2] = setupTableMapping(tbl1, tbl2, crossfields, varargin)
%
% DESCRIPTION: First we describe the structure of the indexing table, then the
% construction of the new indexing table tbl which is made from the given
% crossfields and, finally, the mapping map between vectors indexed according to
% index table tbl1 to vectors indexed according to index table tbl2.
%
% 1) Format of an index table: The index table is a matlab structure whose
% fields are vectors of indices, all the vectors beeing of the same dimension. A
% special field "num" gives the size of the vectors. In this way, we obtain
% multidimensional indexing which is labeled using the field names. If needed,
% we can add local (or linear) indexing of the table, which corresponds to 1,
% ... , tbl.num (see function addLocInd). For any local index i, we have an
% element which corresponds to various elements in the different indexing given
% by the field name. Example: A table of cell-face elements (for example to
% index half-transmissibilities). We choose the field names as 'cells' and
% 'faces', then for each local index i, tbl.cells(i) will give the cell index
% and tbl.faces(i) will give the face index.
%   
% 2) Multi-indexing is difficult because - it not treated carefully - it grows
% exponentially. Here, we provide a way to create, from two given tables tbl1
% and tbl2, new tables that only contains 'relevant' indices and not all
% possible combinations. Let us tbl1 be a table with two fields: tbl1.A and
% tbl1.B, and tbl2 with also two fields tbl2.A and tbl2.C. We choose 'A' as
% crossfield to set up the new table. The result of tbl=setupTableMapping(tbl1,
% tbl2, {'A'}) is the indexing table with fields tbl.A, tbl.B and tbl.C given as
% follows. let index1 = {(i, j) such that i = tbl1.A(m) and j = tbl1.B(m) for some m}
% and similarly index2 = {(i, k) such that i = tbl2.A(n) and j = tbl2.C(n) for
% some m}. Then the table tbl if made of all the triplets (i,j,k) such that
% (i,j) belongs to index1 and (i,k) belongs to index2. Formulated
% differently, we have, for each index l, there exists a unique m and n such
% that tbl.A(l) = tbl1.A(m) = tbl2.A(n), tbl.B(l) = tbl1.B(m) and tbl.C(l) = tbl2.C(n).
%
% 3) The map created by setupTableMapping corresponds to dispactching or
% reduction between the two tables tbl1 and tbl2 across the field given by
% crossfields. Following the notation introduced in 2), the mapping map is a
% linear mapping from vectors of dimension N1=tbl1.num to vectors of dimension
% N2=tbl2.num. For any I = {1, ..., N1} there exist (i1, j1) in index1 and,
% similarly, for any J = {1, ...., N2} there exist (i2, k2) in index2. Then
% map(I, J) = 1 if i1 = i2 and zero otherwise. 
%
%
% The input parameter crossfield can be an empty string or a list of fields. If
% it is empty (crossfields = {''}), then we will get a pure cross-product:
% elements are all duplicated and the dimension of the resulting table is
% tbl.num = tbl1.num*tbl2.num. If it is a list of fields (crossfields = {'A1',
% 'A1'}), then in 2) the comparison is done with respect to all the fields ( for
% each index l, there exists a unique m and n such that tbl.A1(l) = tbl1.A1(m) =
% tbl2.A1(n), tbl.A2(l) = tbl1.A2(m) = tbl2.A2(n), tbl.B(l) = tbl1.B(m) and
% tbl.C(l) = tbl2.C(n).)
%
% PARAMETERS:
%   tbl1        - Indexing table
%   tbl2        - Indexing table
%   crossfields - list of field names used in the way described above. In general it is a list of fields
%                 {'A1', 'A2', 'A3'}. It can be empty {''}. On of the field can be a pair
%                 of fields as for example {'A1', {'A2', 'B'}, 'C'}. Then the field 'A2' of tbl1
%                 and 'B' of tbl2 are considered as the same name. In this case, the first
%                 table tbl1 takes precedence in the naming, that is, the result table tbl
%                 will have field tbl.A2 (and no field tbl.B).
%
%   varargin    - 'crossextend' : This option can be used when  fields in
%   tbl1 and tbl2 have the same name but we want to consider them as
%   different. The syntax for this option is {'A', {'A1', 'A2'}} where 'A' is
%   the common field in tbl1 and tbl2. Then, internally, we replace tbl1.A
%   with tbl1.A1 and tbl2.A with tbl2.A2 and proceed as in the standard case.
%
% KEYWORD ARGUMENTS:
%   
% RETURNS:
%   map  - map between vectors indexed according to tbl1 to vectors indexed
%   according to tbl
%   tbl  - table constructed from tbl1, tbl2 and the given fields
%   crossfields, as described above
%
% SEE ALSO: `projTable`, `sortTable`, `convertArrayToTable`,
% `convertTableToArray`, `addLocInd`, `replacefield`
%


% TODO: not supported situation, the same fieldname is in "other fields" and
% appears as a replacement field name (see below)
    opt = struct('crossextend', [], 'fastunstable', false);
    opt = merge_options(opt, varargin{:});
    
    if ~isempty(opt.crossextend)
        dups = opt.crossextend;
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
            ind1 = tbl1.(fieldname1);
            ind2 = tbl2.(fieldname2);
            nind1 = numel(ind1);
            nind2 = numel(ind2);
            ind = [ind1; ind2];
            [c, ia, ic]= unique(ind);
            ind1 = ic(1 : nind1);
            ind2 = ic((nind1 + 1) : (nind1 + nind2));
            inds1{ifield} = uint64(ind1); 
            inds2{ifield} = uint64(ind2);
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
            tbl.(fieldname) = tbl1.(fieldname)(ind1, 1);
        end
        
        for ifield = 1 : numel(ofields1)
            fieldname = ofields1{ifield};
            tbl.(fieldname) = tbl1.(fieldname)(ind1, 1);
        end

        for ifield = 1 : numel(ofields2)
            fieldname = ofields2{ifield};
            tbl.(fieldname) = tbl2.(fieldname)(ind2, 1);
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
    if any(strcmp(ofields, 'num'));
        ofields = ofields(~strcmp(ofields, 'num'));
    end    
    ofields = ofields(~ismember(ofields, fields));
end
