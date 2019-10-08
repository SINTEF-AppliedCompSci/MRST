classdef Tensor
   
  properties
   nzvals % vector of the nonzero values
   tbl    % struct with named indices ('table')
  end
  
  methods
     function t = Tensor(nzv, ind_tbl)
        assert(isvector(nzv)) % should not be a matrix!
        t.nzvals = nzv(:);
        if nargin < 2
           % if no index set is provided, create a default one
           t.tbl = struct('ix', (1:numel(nzv))', 'num', numel(nzv));
           return 
        elseif iscell(ind_tbl)
           % make a non-sparse, cartesian index table
           indices = cellfun(@(x) 1:x, ind_tbl(2:2:end), 'UniformOutput', false);
           isets = cell(1, numel(indices));
           [isets{:}] = ndgrid(indices{:});
           t.tbl = struct();
           for i = 1:numel(isets)
              t.tbl.(ind_tbl{2*i-1}) = isets{i}(:);
           end
           t.tbl.num = prod([ind_tbl{2:2:end}]);
        else
           assert(isstruct(ind_tbl))
           t.tbl = ind_tbl;
           fnames = Tensor.index_fields(ind_tbl);
           for fn = fnames(:)'
              t.tbl.(fn{:}) = reshape(t.tbl.(fn{:}), [], 1);
           end
           t.tbl.num = numel(t.tbl.(fnames{1}));
        end
     end
     
     function self = contractIn(self, contract_ixs)
        if ~iscell(contract_ixs)
           contract_ixs = {contract_ixs};
        end
        keep_ixs = setdiff(Tensor.index_fields(self.tbl), contract_ixs);
        tbl = projTable(self.tbl, keep_ixs);
        map = setupTableMapping(self.tbl, tbl, Tensor.index_fields(tbl));
        self.tbl = tbl;
        self.nzvals = map * self.nzvals;
     end
     
     
     function res = combineWith(self, other, contract_along)

        if nargin < 3
           contract_along = ...
               intersect(Tensor.index_fields(self.tbl), ...
                         Tensor.index_fields(other.tbl));
           if ~isempty(contract_along)
              contract_along = {contract_along};
           end
        end
        
       % setting identifying names of index fields to contract
       bname = 'contract_me___'; % hopefully an unique name
       fields1 = {}; fields2 = {}; tmpnames = {};
       for i = 1:numel(contract_along)
          tmpnames{i} = [bname, num2str(i)]; %#ok
          fields1 = [fields1, {contract_along{i}{1}, tmpnames{i}}]; %#ok
          fields2 = [fields2, {contract_along{i}{end}, tmpnames{i}}]; %#ok
       end
       
       if ~isempty(fields1)
          tbl_1 = replacefield(self.tbl, fields1);
       else
          tbl_1 = self.tbl;
       end
       
       if ~isempty(fields2)
          tbl_2 = replacefield(other.tbl, fields2);
       else
          tbl_2 = other.tbl;
       end

       fds = ...
           {setdiff(Tensor.index_fields(tbl_1), tmpnames),
            setdiff(Tensor.index_fields(tbl_2), tmpnames),
            tmpnames};
       
       [res_vals, res_tbl] = ...
           contractTable({self.nzvals, tbl_1}, ...
                         {other.nzvals, tbl_2}, ...
                         fds);
       res = Tensor(res_vals, res_tbl);
    end
    
    function self = changeIndexName(self, oldname, newname)
       cur_fields = fields(self.tbl);
       if ~any(strcmp(cur_fields, oldname))
          error('No index with that name.')
       end
       if any(strcmp(cur_fields, newname))
          error('Already an index with the proposed new name.')
       end
       self.tbl.(newname) = self.tbl.(oldname);
       self.tbl = rmfield(self.tbl, oldname);
    end

    function self = add(self, other)
    % check for compatibility
       if ~Tensor.compatible_tables(self, other)
          error('Tried to add incompatible tensors')
       end
       ifields = Tensor.index_fields(self.tbl);
       self = self.sortIndices(ifields);
       other_sorted = other.sortIndices(ifields);
       
       self.nzvals = self.nzvals + other_sorted.nzvals;
    end
    
    function t = mtimes(self, other)
       t = self.combineWith(other);
    end
    
    function t = plus(self, other)
       t = self.add(other);
    end
    
    function self = sortIndices(self, ixset_order)
       num = self.tbl.num;
       [tbl, I] = Tensor.sort_table(self.tbl, ixset_order);
       tbl.num = num;
       self.nzvals = self.nzvals(I);
       self.tbl = tbl;
    end
    
    function M = asSparseMatrix(self, dims)
       if numel(Tensor.index_fields(self.tbl)) > 2
          error('Cannot show a higher dimensional tensor as a matrix.')
       elseif ~Tensor.is_permutation(dims, Tensor.index_fields(self.tbl))
          error('Input ''dims'' must be a permutation of the index sets')
       end
       if numel(dims) == 1
          M = sparse(max(self.tbl.(dims{1})), 1)
          M(self.tbl.(dims{1})) = self.nzvals;
          return
       end
       assert(numel(dims) == 2)
       M = sparse(self.tbl.(dims{1}), self.tbl.(dims{2}), self.nzvals);
    end
  end  
  
  methods(Static)
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
           if any(t1_tbl_sort.(f{:}) ~= t2_tbl_sort.(f{:}))
              return % iscompat = false
           end
        end
        % all tests passed, index sets are equal
        iscompat = true;
     end
  end
  
end
