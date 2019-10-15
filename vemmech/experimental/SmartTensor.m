classdef SmartTensor
   
   properties
      % Cell array of components.  Each component should be a struct with 
      % the fields: 
      % - coefs (vector of coefs)
      % - ixs (matrix where each row constitute a multiindex)
      % - indexnames (names of indices, one for each row of 'ixs'
      components 
   end
   
   methods
      
      function self = SmartTensor(varargin)
         switch nargin
           case 1
             % input is directly in the form of a component
             assert(isstruct(varargin{1}));
             self.components = {varargin{1}};
           case 2
             % input is in form of a matrix (or vector) and a list of (one or
             % two) index names
             coefs = varargin{1};
             indexnames = varargin{2}';
             self.components = SmartTensor.make_matrix_tensor(coefs, indexnames);
           otherwise
             error('Unsupported arguments to constructor.\n');
         end
         
      end
      
      function self = product(self, other)
         self.components = {self.components{:}, other.components{:}};
      end
      
      function self = contract(self, indices)
         error('contract method not yet implemented');
      end

      function self = expandall(self)
         if numel(self.components) == 1
            return
         end
         error('expandall currently unimplemented for multicomponent case.\n');
         % first, we carry out all the contractions in the most
         % computationally friendly order
         
         % then we compute the actual tensor product of the components
      end
      
      function self = plus(self, other)
         self = SmartTensor.apply_binary_operator(self, other, @plus);
      end
      
      function self = minus(self, other)
         self = SmartTensor.apply_binary_operator(self, other, @minus);
      end
      
      function self = rdivide(self, other)
         self = SmartTensor.apply_binary_operator(self, other, @rdivide);
      end
      
      function self = times(self, other)
         self = SmartTensor.apply_binary_operator(self, other, @times);
      end
      
      function t = mtimes(self, other)
         t = self.product(other);
      end
      
      function self = changeIndexName(self, oldnames, newnames)
         if ~iscell(oldnames)
            oldnames = {oldnames}; newnames = {newnames};
         end
         ixnames = self.indexNames(); % verify that it exist before changing it

         for i = 1:numel(oldnames)
            cur_oldname = oldnames{i};
            cur_newname = newnames{i};
            
            for c = numel(self.components):-1:1
               % search for last appearance of name (previous appearances
               % should cancel out pairwise)
               found = strcmp(cur_oldname, self.components{c}.indexnames);
               if any(found)
                  self.components{c}.indexnames{find(found)} = cur_newname;
                  break;
               end
            end
            if ~any(found)
               error(['Could not change index name because it was not found.\n']);
            end
         end
      end
      
      function M = asMatrix(self, ixnames)
         SPARSE_THRESHOLD = 200;
         self = expandall(self);
         % put input variable 'ixnames' on "standard" form
         if ischar(ixnames)
            ixnames = {{ixnames}};
         else
            assert(iscell(ixnames) && numel(ixnames) <= 2)
            for i = 1:numel(ixnames)
               if ~iscell(ixnames{i})
                  ixnames{i} = {ixnames{i}};
               end
            end
         end
         % check that all index names are covered, and compute the
         % permutation
         perm = SmartTensor.get_permutation(horzcat(ixnames{:}), ...
                                            self.indexNames());
         if numel(ixnames) == 1
            ix = compute_1D_index(self.components{1}.ixs(:, perm));
            M = zeros(max(ix), 1);
            M(ix) = full(self.components{1}.coefs);
         elseif numel(ixnames) == 2
            nix1 = numel(ixnames{1});
            ix1 = SmartTensor.compute_1D_index(self.components{1}.ixs(:, perm(1:nix1)));
            ix2 = SmartTensor.compute_1D_index(self.components{1}.ixs(:, perm(nix1+1:end)));
            M = sparse(ix1, ix2, self.components{1}.coefs);
            if prod(size(M)) < SPARSE_THRESHOLD
               M = full(M); % convenient, for small matrices
            end
         end
      end
      
      function self = toInd(self)
         self = self.expandall()
         self.components{1}.coefs = self.components{1}.coefs * 0 + 1;
      end
      
      function ixnames = indexNames(self)
         ixnames = {};
         for c = self.components
            c = c{:};
            for name = c.indexnames
               name = name{:};
               found = strcmp(name, ixnames);
               if any(found)
                  ixnames = ixnames(~found);
               else
                  ixnames = {ixnames{:}, name};
               end
            end
         end
      end
   
   end % end methods

   methods(Static)
      
      function [t1, t2] = make_tensors_compatible(t1, t2)

         % verify that they have the same index sets
         if ~SmartTensor.is_permutation(t1.indexNames(), t2.indexNames())
            error(['Tensors cannot be made compatible as they have different ' ...
                   'index sets.\n']);
         end
         
         % expand both tables
         t1 = t1.expandall();
         t2 = t2.expandall();

         % make sure the order of indices is the same in both tensors
         perm = SmartTensor.get_permutation(t1.indexNames(), t2.indexNames());
         t2.components{1}.indexnames = t1.components{1}.indexnames;
         t2.components{1}.ixs = t2.components{1}.ixs(:, perm);
         
         % fill in missing indices
         [~, I1] = setdiff(t1.components{1}.ixs, t2.components{1}.ixs);
         [~, I2] = setdiff(t2.components{1}.ixs, t1.components{1}.ixs);
         
         t2.components{1}.ixs = [t2.components{1}.ixs; ...
                                 t1.components{1}.ixs(I1,:)];
         t1.components{1}.ixs = [t1.components{1}.ixs; ...
                                 t2.components{1}.ixs(I2,:)];
         t1.components{1}.coefs = [t1.components{1}.coefs; zeros(size(I2))];
         t2.components{1}.coefs = [t2.components{1}.coefs; zeros(size(I1))];
         
         % sort all indices
         t1.components{1} = SmartTensor.sort_indices(t1.components{1});
         t2.components{1} = SmartTensor.sort_indices(t2.components{1});
         
         % tensors should now have exactly the same indices and thus be compatible
      end
      
      function component = sort_indices(component)
         [component.ixs, I] = sortrows(component.ixs);
         component.coefs = component.coefs(I);
      end
      
      function tensor = apply_binary_operator(t1, t2, op)
         [t1, t2] = SmartTensor.make_tensors_compatible(t1, t2);
         assert(numel(t1.components) == 1); % should also be the case for t2 by now
         
         res = t1.components{1};
         res.coefs = op(t1.components{1}.coefs, t2.components{1}.coefs);
         tensor = SmartTensor(res);
      end
      
      function ix = compute_1D_index(multiix)
         if size(multiix, 2) == 1
            ix = multiix;
            return
         end
         tmp = mat2cell(multiix, size(multiix, 1), ones(1, size(multiix, 2)));
         ix = sub2ind(max(multiix), tmp{:});
      end

      function isperm = is_permutation(cells1, cells2)
         isperm = ...
             (numel(cells1) == numel(cells2)) && ...
             isempty(setdiff(cells1, cells2)) && ...
             isempty(setdiff(cells2, cells1));
      end
      
      function perm = get_permutation(cells1, cells2)
         % check that this actually is a permutation
         if ~SmartTensor.is_permutation(cells1, cells2);
            error('provided indices not a permutation of tensor indices.\n');
         end
         perm = zeros(1, numel(cells1));
         for ix = 1:numel(cells1)
            perm(strcmp(cells1{ix}, cells2)) = ix;
         end
      end

      function components = make_matrix_tensor(coefs, indexnames)

         % if user did not wrap the index name in a cell, do it here
         if ~iscell(indexnames)
            indexnames = {indexnames};
         end
         
         if isscalar(coefs) 
            component.coefs = coefs;
            component.ixs = [];
            assert(numel(indexnames) == 0);
         elseif isvector(coefs)
            component.coefs = coefs;
            component.ixs = (1:numel(coefs))';
            assert(numel(indexnames) == 1);
         else
            assert(ismatrix(coefs));
            component.coefs = coefs(:);
            [i, j] = ind2sub(size(coefs), find(coefs));
            component.ixs = [i, j];
         end
         component.indexnames = indexnames(:)';         
         components = {component};
      end
   end
   
end % end classdef
      
   
