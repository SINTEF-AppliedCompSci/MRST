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

      %@OK
      function self = SmartTensor(varargin)
         switch nargin
           case 1
             % input is directly in the form of a component
             assert(isstruct(varargin{1}));
             self.components = {varargin{1}}; %#ok
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

      %@OK
      function self = product(self, other, only_semiproduct)
         if nargin < 3
            only_semiproduct = false;
         end
         
         if ~only_semiproduct
            % we will mark duplicated indices for contraction
            
            % find common index names, and give them new values (since they should
            % be considered to be summed)
            common_names = intersect(self.indexNames(), other.indexNames());
         
            % give new names to the common indices, to mark them for summation
            for nix = 1:numel(common_names)
               % choose an index name that does not overlap with any already
               % used in either of the involved tensors
               flagged_name1 = self.next_unused_contr_ix_name();
               flagged_name2 = other.next_unused_contr_ix_name();
               % in case the two suggested names are different, make sure we
               % pick the 'largest' one.
               choice = {flagged_name1, flagged_name2};
               [~, res] = sort(choice);
               flagged_name = choice{res(end)};
               
               self = self.changeIndexName(common_names{nix}, flagged_name);
               other = other.changeIndexName(common_names{nix}, flagged_name);
            end
         end
         
         self.components = [self.components, other.components];
      end
      
      function self = contract(self, ixname1, ixname2, only_semiproduct)
      % can be called with one index (contraction) or two indices
      % (contraction or semicontraction)
         if nargin < 3
            % Simple contraction in one index.  Check if index is already
            % marked for semi-contraction.  If so, mark it for contraction,
            % otherwise, contract it directly.
            if self.isSemicontracted(ixname1)
               % mark this index for complete contraction
               newname = self.next_unused_contr_ix_name();
               self = self.changeIndexName(ixname1, newname);
            else
               % a simple contraction in one index.  Go ahead and do it
               % directly.
               comp_ix = self.component_with_ix(ixname1);
               self.components{comp_ix} = ...
                   SmartTensor.contract_one_index(self.components{comp_ix}, ...
                                                  ixname1);
               return
            end
         end

         if nargin < 4
            % unless requested otherwise, do a full contraction
            only_semicontract = false; 
         end            
         
         if only_semicontract
            % setting both indexnames equal flags them for semicontraction
            self = self.changeIndexName(ixname2, ixname1);
         else % full contration
            newname = self.next_unused_contr_ix_name();
            self = self.changeIndexName(ixname1, newname);
            self = self.changeIndexName(ixname2, newname);
         end
         
         % check for immediate (semi)contractions (contractions involving
         % only a single component), and carry them out immediately
         self = self.complete_simple_contractions();
         
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
         self = self.product(other, true);
      end
      
      function t = mtimes(self, other)
         t = self.product(other, false);
      end

      %@OK
      function self = sortIndices(self, ixset_order)
      
         self = self.expandall();

         assert(SmartTensor.is_permutation(ixset_order, self.indexNames()));

         perm = SmartTensor.get_permutation(self.indexNames(), ixset_order);

         [self.components{1}.ixs, I] = sortrows(self.components{1}.ixs, perm);

         self.components{1}.coefs = self.components{1}.coefs(I);
         
      end
      
      %@OK
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
               end
            end
         end
      end

      %@OK
      function self = toInd(self)
         self = self.expandall();
         self.components{1}.coefs = self.components{1}.coefs * 0 + 1;
      end
      
      %@OK
      function [ixnames, cixnames, scixnames] = indexNames(self)
      % return tensor index names.  Contracting and semi-contracting index
      % names are returned as second and third return value.  Note that
      % semi-contracting indices are also valid tensor indices.
         ixnames = {};
         cixnames = {};
         scixnames = {};
         for c = self.components
            c = c{:};  %#ok
            for name = c.indexnames
               name = name{:}; %#ok
               if SmartTensor.is_contracting_ix(name)
                  cixnames = [cixnames, name]; %#ok
               else
                  % regular index name
                  ixnames = [ixnames, name]; %#ok
               end
            end
         end
         % determine semi-contracting indices. These are the indices repeated
         % more than once
         [ixnames, ~, ic] = unique(ixnames);
         scixnames = ixnames(find(accumarray(ic, 1) > 1));
         cixnames = unique(cixnames);
      end
      
      %@OK
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
                  ixnames{i} = {ixnames{i}};
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

      function self = expandall(self, expand_tensor)
         if nargin < 2
            expand_tensor = true; % default is true
         end
         
         if numel(self.components) == 1
            return
         end
         
         % get indices flagged for contraction or semi-contraction
         [~, cixs, scixs] = self.indexNames();
         
         ixs = [cixs, scixs];
         while ~isempty(ixs)
            
            % identidy the contraction that requires the least effort
            costs = cellfun(@(x) self.contraction_cost_estimate(x), ixs, ...
                            'uniformoutput', true);
            [~, ix] = find(min(costs)); 
            
            % carry out contraction or semi-contraction
            self = self.execute_semicontraction(ixs{ix});
            
            if any(strcmp(ixs{ix}, cixs))
               % this is a full contraction
               self = self.contract(ixs{ix});
            else 
               % this was a semi-contraction.  We should restore the original
               % index name
               orig_name = split(ixs{ix}, '|');
               assert(numel(orig_name) == 2);
               orig_name = orig_name(2);
               self = self.changeIndexName(ixs{ix}, orig_name);
            end

            % remove processed entry
            ixs = ixs(~strcmp(ixs{ix}, ixs));
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
      
      
      % ---- The following methods represent implementation details, and are not ----
      % -------------------- intended to be directly run by user --------------------

      %@OK
      function self = complete_simple_contractions()
         
         for cix = 1:numel(self.components)
            comp = self.components{cix};
            % check for index names that occur more than once 
            [repnames, ~, ic] = unique(comp.indexnames);
            repnames = repnames(find(accumarray(ic, 1) > 1));
            
            % check if these names occur in any other component; if so, we
            % leave them be, otherwise we contract or semicontract them
            for rname = repnames 
               rname = rname{:};
               found_elsewhere = false;
               for cix2 = 1:numel(self.components)
                  if (cix2 ~= cix) && ...
                     any(strcmp(rname, self.components{cix2}.indexnames))
                     found_elsewhere = true;
                  end
               end
               if ~found_elsewhere 
                  % semi-contract
                  comp = SmartTensor.semi_contract_index(comp, rname);
                  
                  if SmartTensor.is_contracting_ix(rname)
                     % complete the contraction
                     comp = SmartTensor.contract_one_index(comp, rname);
                  end
               end
            end
            self.components{cix} = comp;
         end
      end
      
      function self = execute_semicontraction(self, ixname)
         comp_ixs = self.component_with_ix(ixname, true);
         
         if numel(comp_ixs) == 1
            % the (semi)contraction index was found twice on the same
            % component, perhaps as a result of earlier contractions in other
            % indices.  
            self.components{comp_ixs} = ...
                SmartTensor.semi_contract_two_indices(self.components{comp_ixs}, ...
                                                      ixname, ixname);
         else
            assert(numel(comp_ixs) == 2);
            self.components{comp_ixs(1)} = ...
                SmartTensor.semi_contract_two_indices_two_components( ...
                    self.components{comp_ixs(1)}, ...
                    self.components{comp_ixs(2)}, ...
                    ixname, ixname);
            % remove second component from array of components 
            keep_ix = true(1, numel(self.components));
            keep_ix(comp_ixs(2)) = false;
            self.components = self.components(keep_ix);
         end
      end
      
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
      
      function ixname = next_unused_contr_ix_name(self)
         
         basename = SmartTensor.contracting_name_base();
         [~, current_contr_names] = self.indexNames();
         
         count = 1;
         
         while true
            % search until we find an unused name
            ixname = [ basename, num2str(count)];
            if ~any(strcmp(ixname, current_contr_names))
               % we found a unique name
               return
            else
               % keep searching with a higher count value
               count = count + 1;
            end
         end
      end
      
      function comp_ix = component_with_ix(self, ixname, allow_multiple)
         if nargin < 3
            allow_multiple = false; % the usual case
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
      
      function ixname = contracting_name_base()
         ixname = 'contracting_ix__';
      end
      
      function yesno = is_contracting_ix(ixname)
         % returns true if ixname starts with the contracting name root
         cname = SmartTensor.contracting_name_base();
         yesno = (numel(cname) < numel(ixname)) && ...
                 strcmp(cname, ixname(1:numel(cname)));
      end
      
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

      %@OK
      function comp = semi_contract_index(comp, ixname) 
         
         local_ind = strcmp(ixname, comp.indexnames);
         ixrows = comp.ixs(:, local_ind);
         % keep all entries where the semi-contracted indices have the same value
         keep_ind = sum(abs(diff(ixrows, [], 2)), 2) == 0;
         
         keep_ixvals = ixrows(keep_ind, 1);
         comp.ixs = comp.ixs(keep_ind,:);
         comp.ixs = comp.ixs(:, ~local_ind);
         comp.ixs = [comp.ixs, keep_ixvals];
         comp.coefs = comp.coefs(keep_ind);
         comp.indexnames = [comp.indexnames(~local_ind), ixname];
         
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
