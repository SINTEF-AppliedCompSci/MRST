function [tbl, indstruct] = crossIndexArray(tbl1, tbl2, crossfields, varargin)
%
%
% SYNOPSIS:
%   function [tbl, indstruct] = crossIndexArray(tbl1, tbl2, crossfields, varargin)
%
% DESCRIPTION: 
%   
% Create an IndexArray from two given IndexArrays as the cross-product. The
% result depends on the index that we want to merge which are given in the cell
% array crossfields.
%     
% Multi-indexing is difficult because - if not treated carefully - it grows
% exponentially. Here, we provide a way to create, from two given IndexArrays
% tbl1 and tbl2, new tables that only contains 'relevant' indices and not all
% possible combinations. Let tbl1 be a IndexArray with two fields: tbl1.A and
% tbl1.B, and tbl2 with also two fields tbl2.A and tbl2.C. We choose 'A' as
% crossfield to set up the new table. The result of tbl = crossIndexArray(tbl1,
% tbl2, {'A'}) is the indexing table with fields tbl.A, tbl.B and tbl.C given as
% follows. let index1 = {(i, j) such that i = tbl1.A(m) and j = tbl1.B(m) for
% some m} and similarly index2 = {(i, k) such that i = tbl2.A(n) and j =
% tbl2.C(n) for some m}. Then the IndexArray tbl is made of all the triplets
% (i,j,k) such that (i,j) belongs to index1 and (i,k) belongs to
% index2. Formulated differently, we have that , for each index l, there exists
% a unique m and n such that tbl.A(l) = tbl1.A(m) = tbl2.A(n), tbl.B(l) =
% tbl1.B(m) and tbl.C(l) = tbl2.C(n).
%
% The output indstruct is used to set up mappings between the IndexArray tbl1,
% tbl2 and tbl.
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
%   tbl1        - IndexArray
%   tbl2        - IndexArray
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
%   tbl  - IndexArray constructed from tbl1, tbl2 and the given fields
%   crossfields, as described above 
%   indstruct - structure used to setup mappings between tbl1, tbl2 and tbl3 
%
% SEE ALSO: `IndexArray`
%   `projIndexArray`, `sortIndexArray`, `addLocInd`, `replacefield`.

%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}

% TODO: not supported situation, the same fieldname is in "other fields"
% and appears as a replacement field name (see below)
    opt = struct('crossextend'   , []   , ...
                 'optpureproduct', false, ...
                 'virtual'       , false);
    
    opt = merge_options(opt, varargin{:});
    
    if ~isempty(opt.crossextend)
        dups = opt.crossextend;
        ndups = numel(dups);
        fdnames1 = tbl1.fdnames;
        fdnames2 = tbl2.fdnames;
        for i = 1 : ndups
            dup = dups{i};
            
            fd = dup{1};
            fd1 = dup{2}{1};            
            fd2 = dup{2}{2};

            [isfound1, fdind1] = ismember(fd, fdnames1);
            [isfound2, fdind2] = ismember(fd, fdnames2);
            
            assert(isfound1&isfound2, 'field name not found');
            
            fdnames1{fdind1} = fd1;
            fdnames2{fdind2} = fd2;
        end
        tbl1.fdnames = fdnames1;
        tbl2.fdnames = fdnames2;
    end

    fdnames1 = tbl1.fdnames;
    fdnames2 = tbl2.fdnames;
    
    cfdnames = crossfields; %alias
    nfds = numel(cfdnames);
    
    cfdnames1 = cell(nfds, 1);
    cfdnames2 = cell(nfds, 1);
    
    for ifield = 1 : nfds
        fd = cfdnames{ifield};
        if iscell(fd)
            % handle case with change in replacement field name.
            cfdnames1{ifield} = fd{1};
            cfdnames2{ifield} = fd{2};
        else
            cfdnames1{ifield} = fd;            
            cfdnames2{ifield} = fd;                        
        end
    end
    
    
    n1 = numel(fdnames1);
    [lia, ctblinds1] = ismember(cfdnames1, fdnames1);
    assert(all(lia), 'fieldname not found.');
    otblinds1 = true(n1, 1);
    otblinds1(ctblinds1) = false;
    otblinds1 = find(otblinds1);
    
    n2 = numel(fdnames2);
    [lia, ctblinds2] = ismember(cfdnames2, fdnames2);
    assert(all(lia), 'fieldname not found.');
    otblinds2 = true(n2, 1);
    otblinds2(ctblinds2) = false;
    otblinds2 = find(otblinds2);
    
    
    cfdnames1 = fdnames1(ctblinds1);
    ofdnames1 = fdnames1(otblinds1);
    cfdnames2 = fdnames2(ctblinds2);
    ofdnames2 = fdnames2(otblinds2);    

    [~, otblinds2] = ismember(ofdnames2, fdnames2);
    
    allfdnames = {fdnames1{:}, ofdnames2{:}};
    
    parents.tbl1 = tbl1;
    parents.tbl2 = tbl2;
    
    no2 = numel(ofdnames2);
    n   = numel(allfdnames);
    
    parents.tblinds1 = [(1 : n1)'; zeros(no2, 1)];
    
    tblinds2               = zeros(n, 1);
    tblinds2(ctblinds1)    = ctblinds2;
    tblinds2((n1 + 1) : n) = 1;
    parents.tblinds2 = tblinds2;
    
    n1 = tbl1.num; 
    n2 = tbl2.num;
    
    if (opt.optpureproduct | opt.virtual)
        assert(isempty(crossfields), 'cannot use optpureproduct option in this case');
        % check fieldnames do not intersect
        isnotok = any(ismember(fdnames1, fdnames2)) | any(ismember(fdnames2, ...
                                                          fdnames1));
        assert(~isnotok, 'field names are common, cannot use optpureproduct option in this case');
        
        fdnames = {fdnames1{:}, fdnames2{:}};
        n1 = tbl1.num;
        n2 = tbl2.num;

        if opt.virtual
            tbl = IndexArray([], 'fdnames', fdnames, 'isvirtual', true, 'num', n1*n2);
            return
        end
        
        mat2 = repmat(tbl2.inds, n1, 1);
        mat1 = rldecode(tbl1.inds, n2*ones(n1, 1));
        mat = [mat1, mat2];
        
        tbl = IndexArray([], 'fdnames', fdnames, 'inds', mat);
        
        if nargout > 1
            error(['output not implemented for pureproduct option. It should ' ...
                   'not be difficult though']);
        end
        
        return
        
    end
    mat1 = tbl1.inds(:, ctblinds1);
    mat2 = tbl2.inds(:, ctblinds2);
    
    linind1 = (1 : n1)';
    linind2 = (1 : n2)';
    linind = [linind1; linind2];
    
    mat = [[mat1(:, 1 : nfds), 1*ones(n1, 1)]; ...
           [mat2(:, 1 : nfds), 2*ones(n2, 1)]];
    
    [w, I] = sortrows(mat);
    linind = linind(I);
    
    q = diff(w);
    
    qb = all(q(:, 1 : (end - 1)) == 0, 2);
    qb = qb & (q(:, end) == 1); % pick up element when table index goes from
                                % 1 to 2.
    qb = find(qb);
    % Now qb gives the index in w pointing the last element of from the first table
    % in the series of identical elements. 
    
    [~, n] = rlencode(w, 1);
    nb = rldecode(n, n);
    % Now nb gives for each element in w the number of times this element is
    % repeated.
    
    nb1 = nb(qb);
    nb2 = nb(qb + 1);
    % Now nb1 and nb2 gives, for a series of identical elements, the number
    % that belong to the first and second table, respectively.
    
    ind1 = mcolon(qb - (nb1 - 1), qb);
    ind2 = mcolon(qb + 1, qb + nb2);
    
    
    linind1 = linind(ind1);
    linind2 = linind(ind2);

    [crossind1, crossind2] = crossindconstruct(nb1, nb2);
    
    linind1 = linind1(crossind1);
    linind2 = linind2(crossind2);    

    crossmat = tbl1.inds(linind1, :);
    crossmat = [crossmat, tbl2.inds(linind2, otblinds2)];

    tbl = IndexArray([]);
    tbl.fdnames = allfdnames;
    tbl.inds = crossmat;
    tbl.parents = parents;
    
    indstruct = cell(2, 1);
    indstruct{1}.inds = linind1;
    indstruct{1}.num  = n1;
    indstruct{2}.inds = linind2;
    indstruct{2}.num  = n2;    
    
end


