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
        else
           t.tbl = ind_tbl;
           fnames = fields(ind_tbl);
           for fn = fnames(:)'
              t.tbl.(fn{:}) = reshape(t.tbl.(fn{:}), [], 1);
           end
           t.tbl.num = numel(t.tbl.(fnames{1}));
        end
     end
     
     function res = combineWith(self, other, contract_along)

        if nargin < 3
           contract_along = ...
               setdiff(intersect(fields(self.tbl), fields(other.tbl)), ...
                       {'ind', 'num'});
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
           {setdiff(fields(tbl_1), [tmpnames, {'ind'}, {'num'}]), ...
            setdiff(fields(tbl_2), [tmpnames, {'ind'}, {'num'}]), ...
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
    
    function t = mtimes(self, other)
       t = self.combineWith(other);
    end
    
    function M = asSparseMatrix(self, rowdims, coldims)
    % implement me
       M=sparse(1,1);
    end
  end  
end
