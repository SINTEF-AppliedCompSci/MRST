classdef SmartTensor
   
   properties
      % Cell array of components.  Each component should be a struct with 
      % the fields: 
      % - coefs (vector of coefs)
      % - ixs (matrix where each row constitute a multiindex)
      % - indexnames (names of indices, one for each row of 'ixs')
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
           case 3
             % input is in the form of coefs, ixs and indexnames.  Coefs may
             % be empty, for indicator tensors
             comp.indexnames = varargin{3};
             comp.ixs = varargin{2};
             comp.coefs = varargin{1};
             if isempty(comp.coefs)
                comp.coefs = ones(size(comp.ixs, 1), 1);
             end
             self = SmartTensor(comp);
           otherwise
             error('Unsupported arguments to constructor.');
         end
         
      end
      
      function self = product(self, other)
         
         % find common index names, and give them new values (since they should
         % be considered to be summed)
         common_names = intersect(self.indexNames(), other.indexNames());
         
         % give new names to the common indices, to mark them for summation
         for nix = 1:numel(common_names)
            dummy_name = self.next_unused_dummy_name();
            self = self.changeIndexName(common_names{nix}, dummy_name);
            other = other.changeIndexName(common_names{nix}, dummy_name);
         end
         
         self.components = [self.components, other.components];
      end
      
      function self = semi_contract(self, ixname1, ixname2)
         
         if nargin == 2
            % this should only happen as part of a call to expandall; it is
            % not something a user should do
            comp_ixs = self.component_with_ix(ixname1, true);
            ixname2 = ixname1;
            assert(numel(comp_ixs) == 2);
            comp_ix_1 = comp_ixs(1);
            comp_ix_2 = comp_ixs(2);
         else
            comp_ix_1 = self.component_with_ix(ixname1);
            comp_ix_2 = self.component_with_ix(ixname2);
         end
         
         if comp_ix_1 == comp_ix_2
            self.components{comp_ix_1} = ...
                SmartTensor.semi_contract_two_indices(self.components{comp_ix_1}, ...
                                                      ixname1, ixname2);
         else
            self.components{comp_ix_1} = ...
                SmartTensor.semi_contract_two_indices_two_components( ...
                    self.components{comp_ix_1}, ...
                    self.components{comp_ix_2}, ...
                    ixname1, ixname2);
            % remove comp_ix_2 from array of components (it has been combined
            % with comp_ix_1
            keep_ix = true(1, numel(self.components));
            keep_ix(comp_ix_2) = false;
            self.components = self.components(keep_ix);
         end
         
      end
      
      function self = semi_contract_with(self, other)
         
         common_ixs = intersect(self.indexNames(), other.indexNames());
         common_ixs_dnames = cell(size(common_ixs));
         for i = 1:numel(common_ixs)
            name = common_ixs{i};
            
            dname1 = self.next_unused_dummy_name();
            dname2 = other.next_unused_dummy_name();
            
            if strcmp(dname1, dname2) 
               dname = dname1;
            else
               % choose the highest one
               choice = {dname1, dname2};
               [~, res] = sort(choice);
               dname = choice{res(end)};
            end
            
            self.changeIndexName(name, dname);
            other.changeIndexName(name, dname);
            common_ixs_dnames{i} = dname;
         end
         self = self * other;

         while ~isempty(common_ixs_dnames)
            % identify contraction requiring least effort
            costs = cellfun(@(x) self.contraction_cost_estimate(x), ...
                            common_ixs_dnames, 'uniformoutput', true);
            [~, ix] = find(min(costs));
            
            % carry out semi-contraction
            self = self.semi_contract(common_ixs_dnames{ix});
            self = self.changeIndexName(common_ixs_dnames{ix}, ...
                                        common_ixs{ix});
            
            % remove processed entry
            remove_ind = ~strcmp(common_ixs_dnames{ix}, common_ixs_dnames);
            common_ixs_dnames = common_ixs_dnames(remove_ind);
            common_ixs = common_ixs(remove_ind);
         end
      end
      
      
      function self = contract_one(self, ixname)
         comp_ix = self.component_with_ix(ixname);
         self.components{comp_ix} = ...
             SmartTensor.contract_one_index(self.components{comp_ix}, ixname);
      end
      
      function self = contract_two(self, ixname1, ixname2)
         if nargin == 2
            % usually, the function should be called with two indices.  The
            % present case should only be called from 'expandall', as part of
            % contracting tensors in indices marked for contaction.
            self = self.semi_contract(ixname1);
         else
            self = self.semi_contract(ixname1, ixname2);
         end
         
         self = contract_one(self, ixname1);
      end

      function self = expandall(self, expand_tensor)
         if nargin < 2
            expand_tensor = true; % default is true
         end
         
         if numel(self.components) == 1
            return
         end
         
         cixs = contracting_indices(self); 
         
         while ~isempty(cixs)
            
            % identidy the contraction that requires the least effort
            costs = cellfun(@(x) self.contraction_cost_estimate(x), cixs, ...
                            'uniformoutput', true);
            [~, ix] = find(min(costs)); 
            
            % carry out the contraction
            self = self.contract_two(cixs{ix}); 
                           
            % remove processed entry
            cixs = cixs(~strcmp(cixs{ix}, cixs));
         end
         
         if expand_tensor
            expanded_comp = self.components{1};
            for i = 2:numel(self.components)
               expanded_comp = SmartTensor.tensor_product(expanded_comp, ...
                                                          self.components{i});
               self = SmartTensor(expanded_comp);
            end
         end
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
      
      function self = mpower(self, other)
         self = self.semi_contract_with(other);
      end
      
      function t = mtimes(self, other)
         t = self.product(other);
      end

      function self = sortIndices(self, ixset_order)
      
         self = self.expandall();

         assert(SmartTensor.is_permutation(ixset_order, self.indexNames()));

         perm = SmartTensor.get_permutation(self.indexNames(), ixset_order);

         [self.components{1}.ixs, I] = sortrows(self.components{1}.ixs, perm);

         self.components{1}.coefs = self.components{1}.coefs(I);
         
      end
      
      function self = changeIndexName(self, oldnames, newnames)
         if ~iscell(oldnames)
            oldnames = {oldnames}; newnames = {newnames};
         end

         for i = 1:numel(oldnames)
            cur_oldname = oldnames{i};
            cur_newname = newnames{i};
            
            for c = 1:numel(self.components)
               found = strcmp(cur_oldname, self.components{c}.indexnames);
               if any(found)
                  self.components{c}.indexnames{found} = cur_newname;
                  break;
               end
            end
            if ~any(found)
               error('Could not change index name because it was not found.');
            end
         end
      end
      
      function M = asMatrix(self, ixnames)
         SPARSE_THRESHOLD = 200;
         self = self.expandall();
         
         % if the tensor is an intrinsic scalar, print its value
         if numel(self.indexNames()) == 0
            M = self.components{1}.coefs;
            return
         end
         
         % put input variable 'ixnames' on "standard" form
         if ischar(ixnames)
            ixnames = {{ixnames}};
         else
            assert(iscell(ixnames) && numel(ixnames) <= 2)
            for i = 1:numel(ixnames)
               if ~iscell(ixnames{i})
                  ixnames{i} = {ixnames{i}}; %#ok
               end
            end
         end
         % check that all index names are covered, and compute the
         % permutation
         perm = SmartTensor.get_permutation(self.indexNames(),...
                                            horzcat(ixnames{:}));
                                            
         if numel(ixnames) == 1
            ix = SmartTensor.compute_1D_index(self.components{1}.ixs(:, perm));
            M = zeros(max(ix), 1);
            M(ix) = full(self.components{1}.coefs);
         elseif numel(ixnames) == 2
            nix1 = numel(ixnames{1});
            ix1 = SmartTensor.compute_1D_index(self.components{1}.ixs(:, perm(1:nix1)));
            ix2 = SmartTensor.compute_1D_index(self.components{1}.ixs(:, perm(nix1+1:end)));
            M = sparse(ix1, ix2, self.components{1}.coefs);
            if numel(M) < SPARSE_THRESHOLD
               M = full(M); % convenient, for small matrices
            end
         end
      end
      
      function self = toInd(self)
         self = self.expandall();
         self.components{1}.coefs = self.components{1}.coefs * 0 + 1;
      end
      
      function ixnames = indexNames(self)
         ixnames = {};
         for c = self.components
            c = c{:};  %#ok
            for name = c.indexnames
               name = name{:}; %#ok
               found = strcmp(name, ixnames);
               if any(found)
                  % this procedure ensures cancellation of indice names meant
                  % for summation
                  ixnames = ixnames(~found);
               else
                  ixnames = [ixnames, name]; %#ok
               end
            end
         end
      end
      
      % ---- The following methods represent implementation details, and are not ----
      % -------------------- intended to be directly run by user --------------------

      function cost = contraction_cost_estimate(self, cix)

         comps = self.components;
         comp_ixs = component_with_ix(self, cix, true);
         assert(numel(comps) == 2);
         
         cost = 0;
         for i = comp_ixs
            % we consider the cost to be proportional to the size of the
            % tensor, were it to be fully expanded
            cost = cost + size(comps{i}.ixs, 1) ^ size(comps{i}.ixs, 2);
         end         
         
      end
      
      function cixs = contracting_indices(self)
         regular_indices = indexNames(self);
         
         % contracting indices are those that are not regular indices
         % (i.e. they are repeated twice)
         cixs = {};
         for c = self.components
            c = c{:}; %#ok
            for name = c.indexnames
               name = name{:}; %#ok
               found = strcmp(name, regular_indices);
               if ~any(found)
                  % this is not a regular index.  It must be a contracting
                  % index
                  cixs = [cixs, name]; %#ok
               end
            end
         end
         cixs = unique(cixs);
      end
      
      function dummy_name = next_unused_dummy_name(self)
         
         current_index_names = self.indexNames();
         
         count = 1;
         template = 'summing_ix_%i__';
         while true
            % search until we find an unused dummy name
            dummy_name = sprintf(template, count);
            if ~any(strcmp(dummy_name, current_index_names))
               % we found a unique name
               return
            else
               % keep searching with a higher count value
               count = count + 1;
            end
         end
      end
      
      function check_valid_ix(self, ixname)
         if ~any(strcmp(ixname, self.indexNames()))
            error('Unknown index.');
         end
      end
      
      function comp_ix = component_with_ix(self, ixname, allow_multiple)
         if nargin < 3
            allow_multiple = false; % the usual case
         end
         if ~allow_multiple
            self.check_valid_ix(ixname);
         end

         res = [];
         
         for comp_ix = numel(self.components):-1:1
            if any(strcmp(ixname, self.components{comp_ix}.indexnames))
               if allow_multiple
                  res = [res, comp_ix]; %#ok
               else
                  return;
               end
            end
         end
         comp_ix = res; % multiple compoents (or zero)
      end
      
   end % end methods

   methods(Static)
      
      function tensor = ind(tensor)
         tensor = tensor.toInd();
      end
         
      function comp = tensor_product(comp1, comp2)
         
         if numel(comp1.indexnames) == 0
            % comp1 is an intrinsic scalar
            comp = comp2;
            comp.coefs = comp1.coefs(1) * comp.coefs;
            return 
         elseif numel(comp2.indexnames) == 0
            % comp2 is an intrinsic scalar
            comp = comp1;
            comp.coefs = comp2.coefs(1) * comp.coefs;
            return
         else
            % both components were nontrivial tensors.  Expand tensor product
            comp.indexnames = [comp1.indexnames, comp2.indexnames];
            comp.coefs = kron(comp1.coefs, comp2.coefs);
            comp.ixs = [repelem(comp1.ixs, size(comp2.ixs, 1), 1), ...
                        repmat(comp2.ixs, size(comp1.ixs, 1), 1)];
         end
      end
                                                                        
      function comp = semi_contract_two_indices_two_components(comp1, comp2, ...
                                                           ixname1, ixname2)
      % carry out semi-contraction, involving two components
         
         % eliminate entries that do not contribute
         ix1 = find(strcmp(ixname1, comp1.indexnames));
         ix2 = find(strcmp(ixname2, comp2.indexnames));
         
         common_indices = intersect(comp1.ixs(:, ix1), comp2.ixs(:, ix2));
         num_unique_ix = numel(common_indices);
         keep1 = ismember(comp1.ixs(:, ix1), common_indices);
         keep2 = ismember(comp2.ixs(:, ix2), common_indices);
         
         comp1.ixs = comp1.ixs(keep1, :);  comp1.coefs = comp1.coefs(keep1);
         comp2.ixs = comp2.ixs(keep2, :);  comp2.coefs = comp2.coefs(keep2);
         
         % identify location of each relevant index (after sorting entries)
         [comp1.ixs, I1] = sortrows(comp1.ixs, ix1); comp1.coefs = comp1.coefs(I1);
         [comp2.ixs, I2] = sortrows(comp2.ixs, ix2); comp2.coefs = comp2.coefs(I2);
         
         start1 = [1; find(diff(comp1.ixs(:, ix1))) + 1];
         end1   = [start1(2:end)-1; size(comp1.ixs, 1)];

         start2 = [1; find(diff(comp2.ixs(:, ix2))) + 1];
         end2   = [start2(2:end)-1; size(comp2.ixs, 1)];
         
         % split up matrices
         ixsplit1 = arrayfun(@(x) comp1.ixs(start1(x):end1(x), :), 1:num_unique_ix, ...
                             'uniformoutput', false);
         ixsplit2 = arrayfun(@(x) comp2.ixs(start2(x):end2(x), :), 1:num_unique_ix, ...
                             'uniformoutput', false);
         coefsplit1 = arrayfun(@(x) comp1.coefs(start1(x):end1(x)), 1:num_unique_ix, ...
                               'uniformoutput', false);
         coefsplit2 = arrayfun(@(x) comp2.coefs(start2(x):end2(x)), 1:num_unique_ix, ...
                               'uniformoutput', false);
         
         % function producing all permutation of rows in matrices u and v
         loc_kron = @(u, v) [repmat(u, size(v,1), 1), repelem(v, size(u, 1), 1)];
         
         % make array with all indices that are to be multiplied together
         ixarr = arrayfun(@(x) loc_kron(ixsplit1{x}, ixsplit2{x}), 1:num_unique_ix, ...
                          'uniformoutput', false);
         ixarr = vertcat(ixarr{:});
         
         % produce corresponding permutation also for the coefficient vectors
         % (which may be ADI, so we must be a bit careful)
         
         cfarr1 = arrayfun(@(x) repmat(coefsplit1{x}, numel(coefsplit2{x}),1), ...
                           1:num_unique_ix, 'uniformoutput', false);
         cfarr1 = vertcat(cfarr1{:});
         
         cfarr2 = arrayfun(@(x) coefsplit2{x}(reshape(repmat(1:numel(coefsplit2{x}),...
                                                           numel(coefsplit1{x}), 1), ...
                                              [], 1)), ...
                           1:num_unique_ix, 'uniformoutput', false);
         cfarr2 = vertcat(cfarr2{:});

         % constructing result
         remove_ix = find(strcmp(ixname2, comp2.indexnames)) + numel(comp1.indexnames);
         
         comp.indexnames = [comp1.indexnames, ...
                            comp2.indexnames(~strcmp(ixname2, comp2.indexnames))];
         comp.coefs = cfarr1 .* cfarr2;
         comp.ixs = ixarr;
         comp.ixs(:, remove_ix) = [];

      end
      
      function comp = semi_contract_two_indices(comp, ixname1, ixname2)

         local_ind_1 = strcmp(ixname1, comp.indexnames);
         local_ind_2 = strcmp(ixname2, comp.indexnames);
         
         ixrows = comp.ixs(:, local_ind_1 | local_ind_2);
         assert(size(ixrows, 2) == 2);
         
         keep_ixs = ixrows(:,1) == ixrows(:,2);
         
         comp.ixs = comp.ixs(keep_ixs, :); % remove discarded nonzero entries
         comp.ixs = comp.ixs(:, ~local_ind_2); % remove second index
         
         comp.indexnames = comp.indexnames(~local_ind_2); % remove second indexname
      
         comp.coefs = comp.coefs(keep_ixs);
      end
      
      function comp = contract_one_index(comp, ixname)
         local_ind = strcmp(ixname, comp.indexnames);

         % get rid of contracting index
         comp.indexnames = comp.indexnames(~local_ind);
         comp.ixs = comp.ixs(:, ~local_ind);
         if isempty(comp.ixs)
            % result is an intrinsic scalar
            comp.coefs = sum(comp.coefs);
            comp.ixs = []; % get rid of "ghost" dimensions
            return
         end
         
         % sum up other elements at the correct place
         map = accumarray([SmartTensor.compute_1D_index(comp.ixs), ...
                           (1:size(comp.ixs,1))'], ...
                          1, [], [], [], true);
         comp.coefs = map * comp.coefs;
         
         % extract nonzeros and recompute indices
         nz = find(comp.coefs);
         comp.coefs = comp.coefs(nz);
         
         logical_size = max(comp.ixs);
         reindex = cell(size(comp.ixs, 2), 1);
         
         [reindex{:}] = ind2sub(logical_size, nz);
         
         comp.ixs = [reindex{:}];
      end
      
      function [t1, t2] = make_tensors_compatible(t1, t2)

         % verify that they have the same index sets
         if ~SmartTensor.is_permutation(t1.indexNames(), t2.indexNames())
            error(['Tensors cannot be made compatible as they have different ' ...
                   'index sets.']);
         end
         
         % expand both tables
         t1 = t1.expandall();
         t2 = t2.expandall();

         % make sure the order of indices is the same in both tensors
         perm = SmartTensor.get_permutation(t2.indexNames(), t1.indexNames());
         t2.components{1}.indexnames = t1.components{1}.indexnames;
         t2.components{1}.ixs = t2.components{1}.ixs(:, perm);
         
         % fill in missing indices
         [~, I1] = setdiff(t1.components{1}.ixs, t2.components{1}.ixs, 'rows');
         [~, I2] = setdiff(t2.components{1}.ixs, t1.components{1}.ixs, 'rows');
         
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
         if ~SmartTensor.is_permutation(cells1, cells2)
            error('provided indices not a permutation of tensor indices.');
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
            nz = find(coefs);
            component.coefs = coefs(nz);
            component.ixs = nz(:);
            assert(numel(indexnames) == 1);
         else
            assert(ismatrix(coefs));
            nz = find(coefs); % only keep nonzeros
            %nz = (1:numel(coefs))';
            [i, j] = ind2sub(size(coefs), nz);
            component.coefs = coefs(nz);
            component.ixs = [i, j];
         end
         component.indexnames = indexnames(:)';         
         components = {component};
      end
   end
   
end % end classdef
