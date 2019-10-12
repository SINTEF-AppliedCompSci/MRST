classdef Tensor
   
  properties
   nzvals % vector of the nonzero values
   tbl    % struct with named indices ('table')
  end
  
  methods

     function self = Tensor(varargin)
     % dispatch constructor based on argument
        tt = [];
        switch nargin
          case 1
            a = varargin{1};
            if Tensor.isTable(a)
               tt = Tensor.make_indicator_tensor(a);
            elseif isvector(a)
               tt = Tensor.make_vector_tensor(a);
            elseif ismatrix(a)
               tt = Tensor.make_matrix_tensor(a);
            end
          case 2
            [a, b] = deal(varargin{:});
            if isvector(a)
               if ischar(b)
                  tt = Tensor.make_vector_tensor(a, b);
               elseif Tensor.isLabelDimList(b)
                  tt = Tensor.make_full_tensor(a, b);
               elseif Tensor.isTable(b)
                  tt = Tensor.make_generic_tensor(a, b);
               end
            elseif ismatrix(a) && Tensor.isLabelList(b)
               tt = Tensor.make_matrix_tensor(a, b);
            end
          otherwise
            % do nothing, we'll report the error further down
        end
        if ~isempty(tt)
           self.nzvals = tt.nzvals;
           self.tbl = tt.tbl;
           return
        end
        % we'll only get here if no valid argument combination was recognized
        error('Unrecognized argument list for Tensor constructor.')
     end

     
     function self = project(self, other)
        common_indexsets = Tensor.find_common_indices(self, other);
        
        % @ temporary fix
        self.tbl.num = numel(self.nzvals);
        other.tbl.num = numel(other.nzvals);
        [~, combined_tbl] = ...
            setupTableMapping(self.tbl, other.tbl, common_indexsets, ...
                              'fastunstable', false);
        
        m1 = setupTableMapping(self.tbl, combined_tbl, ...
                               Tensor.index_fields(self.tbl), 'fastunstable', ...
                               false);
        m2 = setupTableMapping(other.tbl, combined_tbl, ...
                               Tensor.index_fields(other.tbl), 'fastunstable', ...
                               false);
        self.nzvals = (m1 * self.nzvals) .* (m2 * other.nzvals);
        self.tbl = combined_tbl;
        
     end

     function self = contract(self, contract_ixs)
        
        if ~iscell(contract_ixs)
           contract_ixs = {contract_ixs};
        end
        
        self.tbl = rmfield(self.tbl, contract_ixs);
        [ix1D, dims] = Tensor.compute_1D_index(self.tbl);
        map = sparse(ix1D, 1:numel(self.nzvals), 1);
        
        % @@the following line has been known to cause memory issues, so it has
        % been replaced with the equivalent (but slightly more complex) code below
        numel(self.nzvals)
        nzvals_sparse = map * self.nzvals;
        
        % @@ this code replaces the line above
        %nzvals_sparse = sum(bsxfun(@times, map, self.nzvals'), 2);
        
        
        [ix,~, self.nzvals] = find(nzvals_sparse);
        subs = cell(numel(dims), 1); 
        [subs{:}] = ind2sub(dims, ix);
        subs = horzcat(subs{:});
        
        fields = Tensor.index_fields(self.tbl);
        for f = 1:numel(fields)
           self.tbl.(fields{f}) = subs(:, f);
        end
     end
     
     function self = formalProduct(self, other)
        common_indexsets = Tensor.find_common_indices(self, other);
        self = self.project(other).contract(common_indexsets);
     end
          
     function self = changeIndexName(self, oldname, newname)
        if ~iscell(oldname)
           oldname = {oldname}; newname = {newname};
        end
        cur_fields = fields(self.tbl);
        
        for i = 1:numel(oldname)
           cur_oldname = oldname{i};
           cur_newname = newname{i};
           
           if ~any(strcmp(cur_fields, cur_oldname))
              error('No index with that name.')
           end
           if any(strcmp(cur_fields, cur_newname))
              error('Already an index with the proposed new name.')
           end
           
           self.tbl.(cur_newname) = self.tbl.(cur_oldname);
           self.tbl = rmfield(self.tbl, cur_oldname);
        end
     end

    function self = plus(self, other)
       self = Tensor.apply_binary_operator(self, other, @plus);
    end

    function self = minus(self, other)
       self = Tensor.apply_binary_operator(self, other, @minus);
    end

    function self = rdivide(self, other)
       self = Tensor.apply_binary_operator(self, other, @rdivide);
    end

    function self = times(self, other)
       self = Tensor.apply_binary_operator(self, other, @times);
    end       
    
    function t = mtimes(self, other)
       t = self.formalProduct(other);
    end
    
    function self = sortIndices(self, ixset_order)
       num = numel(self.nzvals);
       [new_tbl, I] = Tensor.sort_table(self.tbl, ixset_order);
       self.nzvals = self.nzvals(I);
       self.tbl = new_tbl;
    end
    
    function M = asMatrix(self, dims)
       SPARSE_THRESHOLD = 200;
       % put input variable 'dims' on "standard" form
       if ischar(dims)
          dims = {{dims}};
       else
          assert(iscell(dims) && numel(dims) <= 2)
          for i = 1:numel(dims)
             if ~iscell(dims{i})
                dims{i} = {dims{i}};
             end
          end
       end
       % check that all index names are covered
       assert(Tensor.is_permutation(horzcat(dims{:}), ...
              Tensor.index_fields(self.tbl)));
       if numel(dims) == 1
          ix1 = Tensor.compute_1D_index(self.tbl, dims{1});
          M = zeros(max(ix1), 1);
          M(ix1) = self.nzvals;
       elseif numel(dims) == 2
          ix1 = Tensor.compute_1D_index(self.tbl, dims{1});
          ix2 = Tensor.compute_1D_index(self.tbl, dims{2});
          M = sparse(ix1, ix2, self.nzvals);
          if prod(size(M)) < SPARSE_THRESHOLD
             M = full(M); % convenient, for small matrices
          end
       end
    end
    
    function self = lockIndex(self, ixname, ixval)    
       if ~any(strcmp(ixname, Tensor.index_fields(self.tbl)))
          error('Tensor has no index with that name.')
       end
       keep_ixs = self.tbl.(ixname) == ixval;
       self.nzvals = self.nzvals(keep_ixs);
       tbl = struct();
       for field = Tensor.index_fields(self.tbl)'
          tbl.(field{:}) = self.tbl.(field{:})(keep_ixs);
       end
       tbl = rmfield(tbl, ixname);
       self.tbl = tbl;
    end

     function self = toInd(self)
        self.nzvals = ones(numel(self.nzvals), 1);
     end
    
     function self = pruneZeros(self)
        prune_ix = self.nzvals==0;
        self.nzvals = self.nzvals(~prune_ix);
        for field = Tensor.index_fields(self.tbl)'
           self.tbl.(field{:}) = self.tbl.(field{:})(~prune_ix);
        end
     end
     
  end % end regular, public methods
     
  methods(Static)
     
     function [ix, dims] = compute_1D_index(tbl, dims)
        if nargin > 1
           % keep only the requested dimensions
           remove = setdiff(Tensor.index_fields(tbl), dims);
           tbl = rmfield(tbl, remove);
        end
        
        tbl_ix = cellfun(@(f) tbl.(f), Tensor.index_fields(tbl), ...
                          'uniformoutput', false);
        tbl_ix = [tbl_ix{:}];
        dims = max(tbl_ix);
        stride = cumprod([1, dims(1:end-1)]);
        
        ix = sum(bsxfun(@times, tbl_ix-1, stride), 2) + 1;
        % (could also have used sub2ind to do this)
     end
     
     function tensor = ind(tensor)
        tensor = tensor.toInd();
     end
     
     function flds = index_fields(tbl)
        flds = setdiff(fields(tbl), {'num', 'ind'});
     end
     
     function [tbl, I] = sort_table(tbl, fds)
        ifields = Tensor.index_fields(tbl);
        
        % ensure 'fds' is a permutation of 'ifields'
        assert(Tensor.is_permutation(fds, ifields));
        
        A = convertTableToArray(tbl, fds);
        [A, I] = sortrows(A);
        tbl = convertArrayToTable(A, fds);
     end
     
     function isperm = is_permutation(cellarr1, cellarr2)
        isperm = ...
            (numel(cellarr1) == numel(cellarr2)) && ...
            isempty(setdiff(cellarr1, cellarr2)) && ...
            isempty(setdiff(cellarr2, cellarr1));
     end
     
     function iscompat = compatible_tables(t1, t2)
        fields1 = Tensor.index_fields(t1.tbl);
        fields2 = Tensor.index_fields(t2.tbl);
        iscompat = false;
        
        if ~(Tensor.is_permutation(fields1, fields2))
           return % iscompat = false
        end
        
        t1_tbl_sort = Tensor.sort_table(t1.tbl, fields1);
        t2_tbl_sort = Tensor.sort_table(t2.tbl, fields1); 
        
        for f = fields1'
           m1 = t1_tbl_sort.(f{:});
           m2 = t2_tbl_sort.(f{:});
           if numel(m1) ~= numel(m2)
              return % iscompat = false
           elseif any(m1 ~= m2)
              return % iscompat = false
           end
        end
        % all tests passed, index sets are equal
        iscompat = true;
     end
     
     function [t1, t2] = make_tables_compatible(t1, t2)
        fields1 = Tensor.index_fields(t1.tbl);
        fields2 = Tensor.index_fields(t2.tbl);
        
        if Tensor.compatible_tables(t1, t2)
           % already compatible
           return
        elseif ~(Tensor.is_permutation(fields1, fields2))
           error('Index tables cannot be made compatible.')
        end
        % tables are not compatible, but can be made so
        
        t1_entries = cellfun(@(x) t1.tbl.(x)(:), fields1, 'uniformoutput', false);
        t1_entries = [t1_entries{:}];
        t2_entries = cellfun(@(x) t2.tbl.(x)(:), fields1, 'uniformoutput', false);
        t2_entries = [t2_entries{:}];

        t2_missing = setdiff(t1_entries, t2_entries, 'rows');
        t1_missing = setdiff(t2_entries, t1_entries, 'rows');
        
        for i = 1:numel(fields1)
           fld = fields1{i};
           t1.tbl.(fld) = [t1.tbl.(fld); t1_missing(:, i)];
           t2.tbl.(fld) = [t2.tbl.(fld); t2_missing(:, i)];
        end
        t1.nzvals = [t1.nzvals; zeros(size(t1_missing, 1), 1)];
        t2.nzvals = [t2.nzvals; zeros(size(t2_missing, 1), 1)];
        assert(numel(t1.nzvals) == numel(t2.nzvals));
     end
     
     
     function res = apply_binary_operator(t1, t2, op)   
     
        [t1, t2] = Tensor.make_tables_compatible(t1, t2);
     
        ifields = Tensor.index_fields(t1.tbl);
        res = t1.sortIndices(ifields);
        t2_sorted = t2.sortIndices(ifields);
        
        res.nzvals = op(res.nzvals,  t2_sorted.nzvals);
     end
     
     function isetlist = find_common_indices(t1, t2)
        isetlist = intersect(Tensor.index_fields(t1.tbl), ...
                             Tensor.index_fields(t2.tbl));
     end

     function res = isTable(var)
        res = isstruct(var); % @ a more rigorous test could be applied here
     end
     
     function res = isLabelList(var)
        % check that 'var' is a cell array of labels
        res = iscell(var) && all(cellfun(@ischar, var));
     end
     
     function res = isLabelDimList(var)
     % check that 'var' is a cell array of alternating labels and numbers
        res = iscell(var) && ...
              Tensor.isLabelList(var(1:2:end)) && ...
              all(cellfun(@isscalar, var(2:2:end)));
     end

     function t = make_indicator_tensor(tbl)
        flds = Tensor.index_fields(tbl);
        t = Tensor.make_generic_tensor(ones(numel(tbl.(flds{1})), 1), tbl);
     end
     
     function t = make_vector_tensor(vec, label)
        if nargin < 2
           label = 'i';
        end
        t.nzvals = vec(:);
        t.tbl = struct(label, (1:numel(vec))', 'num', numel(vec));
     end
     
     function t = make_matrix_tensor(mat, labels)
        if nargin < 2
           labels = {'row', 'column'};
        end
        ixs = find(mat);
        t.nzvals = mat(ixs);
        [r, c] = ind2sub(size(mat), ixs);
        t.tbl = struct(labels{1}, r, labels{2}, c, 'num', numel(ixs));
     end
     
     function t = make_full_tensor(coefs, label_dims)
     % make a non-sparse, cartesian index table
        t.nzvals = coefs(:);
        indices = cellfun(@(x) 1:x, label_dims(2:2:end), 'UniformOutput', false);
        isets = cell(1, numel(indices));
        [isets{:}] = ndgrid(indices{:});
        t.tbl = struct();
        for i = 1:numel(isets)
           t.tbl.(label_dims{2*i-1}) = isets{i}(:);
        end
     end
     
     function t = make_generic_tensor(coefs, tbl)
        assert(Tensor.isTable(tbl))
        t.nzvals = coefs(:);
        t.tbl = tbl;
        fnames = Tensor.index_fields(tbl);
        for fn = fnames(:)'
           t.tbl.(fn{:}) = reshape(t.tbl.(fn{:}), [], 1);
        end
     end
  end % end static methods
end 



    %  function res = combineWith(self, other, contract_along)

    %     if nargin < 3
    %        contract_along = ...
    %            intersect(Tensor.index_fields(self.tbl), ...
    %                      Tensor.index_fields(other.tbl));
    %        if ~isempty(contract_along)
    %           contract_along = cellfun(@(x) {x, x}, contract_along, ...
    %                                    'UniformOutput', false);
    %        end
    %     end
        
    %    % setting identifying names of index fields to contract
    %    bname = 'contract_me___'; % hopefully an unique name
    %    fields1 = {}; fields2 = {}; tmpnames = {};
    %    for i = 1:numel(contract_along)
    %       tmpnames{i} = [bname, num2str(i)]; %#ok
    %       fields1 = {fields1{:}, {contract_along{i}{1}, tmpnames{i}}}; %#ok
    %       fields2 = {fields2{:}, {contract_along{i}{end}, tmpnames{i}}}; %#ok
    %    end
       
    %    if ~isempty(fields1)
    %       tbl_1 = replacefield(self.tbl, fields1);
    %    else
    %       tbl_1 = self.tbl;
    %    end
       
    %    if ~isempty(fields2)
    %       tbl_2 = replacefield(other.tbl, fields2);
    %    else
    %       tbl_2 = other.tbl;
    %    end

    %    fds = ...
    %        {setdiff(Tensor.index_fields(tbl_1), tmpnames), ...
    %         setdiff(Tensor.index_fields(tbl_2), tmpnames), ...
    %         tmpnames};
       
    %    [res_vals, res_tbl] = ...
    %        contractTable({self.nzvals, tbl_1}, ...
    %                      {other.nzvals, tbl_2}, ...
    %                      fds);
    %    res = Tensor(res_vals, res_tbl);
    % end
    
     % function t = Tensor(nzv, ind_tbl)
     %    if nargin == 1 && isstruct(nzv)
     %       % special case: user wants an indicator tensor, and nzv is really a
     %       % table struct here.  Let all cofficients be ones
     %       flds = Tensor.index_fields(nzv);
     %       t = Tensor(ones(numel(nzv.(flds{1})), 1), nzv);
     %       return 
     %    end
        
     %    assert(isvector(nzv)) % should not be a matrix!
     %    t.nzvals = nzv(:);
     %    if nargin < 2
     %       % if no index set is provided, create a default one
     %       t.tbl = struct('ix', (1:numel(nzv))', 'num', numel(nzv));
     %       return 
     %    elseif ischar(ind_tbl)
     %       t.tbl = struct(ind_tbl, (1:numel(nzv))', 'num', numel(nzv));
     %       return
     %    elseif iscell(ind_tbl)
     %       % make a non-sparse, cartesian index table
     %       indices = cellfun(@(x) 1:x, ind_tbl(2:2:end), 'UniformOutput', false);
     %       isets = cell(1, numel(indices));
     %       [isets{:}] = ndgrid(indices{:});
     %       t.tbl = struct();
     %       for i = 1:numel(isets)
     %          t.tbl.(ind_tbl{2*i-1}) = isets{i}(:);
     %       end
     %       t.tbl.num = prod([ind_tbl{2:2:end}]);
     %    else
     %       assert(isstruct(ind_tbl))
     %       t.tbl = ind_tbl;
     %       fnames = Tensor.index_fields(ind_tbl);
     %       for fn = fnames(:)'
     %          t.tbl.(fn{:}) = reshape(t.tbl.(fn{:}), [], 1);
     %       end
     %       if isempty(fnames)
     %          t.tbl.num = 1;
     %       else
     %          t.tbl.num = numel(t.tbl.(fnames{1}));
     %       end
     %    end
     % end
     
    % function M = asMatrix(self, dims)
    %    SPARSE_THRESHOLD = 200;
    %    if numel(Tensor.index_fields(self.tbl)) > 2
    %       error('Cannot show a higher dimensional tensor as a matrix.')
    %    elseif ~Tensor.is_permutation(dims, Tensor.index_fields(self.tbl))
    %       error('Input ''dims'' must be a permutation of the index sets')
    %    end
    %    if numel(dims) == 1
    %       M = zeros(max(self.tbl.(dims{1})), 1);
    %       % M = sparse(max(self.tbl.(dims{1})), 1);
    %       M(self.tbl.(dims{1})) = self.nzvals;
    %       return
    %    end
    %    assert(numel(dims) == 2)
    %    M = sparse(self.tbl.(dims{1}), self.tbl.(dims{2}), self.nzvals);
    %    if prod(size(M)) < SPARSE_THRESHOLD
    %       M = full(M); % we return a full matrix when reasonable
    %    end
       
    % end
    
     % function self = contract_old(self, contract_ixs)
     %    if ~iscell(contract_ixs)
     %       contract_ixs = {contract_ixs};
     %    end
     %    keep_ixs = setdiff(Tensor.index_fields(self.tbl), contract_ixs);
     %    projtbl = projTable(self.tbl, keep_ixs); 
     %    map = setupTableMapping(self.tbl, projtbl, Tensor.index_fields(projtbl), ...
     %                            'fastunstable', false);
     %    %@@ NB: The following line is wrong!  The number of entries is equal
     %    %to the number of entries in the original, uncontracted tensor, so
     %    %there will be dupplication.  Function unusable, and kept for future
     %    %reference only.
     %    self.tbl = projtbl;
     %    self.nzvals = map * self.nzvals;
     % end
    
 