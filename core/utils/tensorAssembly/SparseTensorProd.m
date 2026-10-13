classdef SparseTensorProd

    properties (SetAccess = private)

        tensor1 % first tensor argument
        tensor2 % second tensor argument

    end

    properties

        reducefds    % Field names for reduction
        mergefds     % Field names for merging
        
        replacefds1  % Possibility to change field names of first table
                     % (before setting up product)
        replacefds2  % Possibility to change field names of second table
                     % (before setting up product)
        replacefds3  % Possibility to change field names of product table (after setting
                     % up product) for the default tbl3 (not applied if tbl3 is
                     % given as input). This option is little used, as it can be
                     % avoid by using only replacefds1 and replacefds2

        reducefds1   % Field names that belongs to tbl1 and not tbl2 that
                     % will be reduced
        reducefds2   % Field names that belongs to tbl2 and not tbl1 that
                     % will be reduced
        
    end
   
    methods
        
        function prod = SparseTensorProd(tensor1, tensor2, varargin)
            
            prod = merge_options(prod, varargin{:}); 

            prod.tensor1 = tensor1;
            prod.tensor2 = tensor2;
            
            if isempty(prod.reducefds)
                prod.reducefds = {};
            end
            
            if isempty(prod.reducefds1)
                prod.reducefds1 = {};
            end
            
            if isempty(prod.reducefds2)
                prod.reducefds2 = {};
            end

            if isempty(prod.mergefds)
                prod.mergefds = {};
            end            
            
        end

        function C = eval(prod)
            
            reducefds  = prod.reducefds;
            mergefds   = prod.mergefds;

            crossfds  = {reducefds{:}, mergefds{:}};

            if ~isempty(prod.replacefds1)
                prod.tensor1.tbl = replacefield(prod.tensor1.tbl, prod.replacefds1);
            end
            if ~isempty(prod.replacefds2)
                prod.tensor2.tbl = replacefield(prod.tensor2.tbl, prod.replacefds2);
            end

            % some sanity checks on the table's field names.
            fds1 = prod.tensor1.tbl.fdnames;
            fds2 = prod.tensor2.tbl.fdnames;

            assert(all(ismember(mergefds, fds1)), ['There exist merge fields that do ' ...
                                'not belong to fields of first table']);
            assert(all(ismember(mergefds, fds2)), ['There exist merge fields that do ' ...
                                'not belong to fields of second table']);
            ofds1 = fds1(~ismember(fds1, crossfds));
            ofds2 = fds2(~ismember(fds2, crossfds));
            assert(all(~ismember(ofds1, ofds2)) & all(~ismember(ofds2, ofds1)), ...
                   ['There exist fields with same name in first and second ' ...
                    'table that are neither merged or reduced.']);

            if ~isempty(prod.reducefds1)
                
                assert(all(ismember(prod.reducefds1, ofds1)), ['There exist pre-reduced fields of tbl1 ' ...
                                                               'that do not belong to ofds1']);
                stm = SparseTensorMap(prod.tensor1);
                stm.toTbl    = projIndexArray(prod.tensor1.tbl, fds1(~ismember(fds1, prod.reducefds1)));
                stm.mergefds = stm.toTbl.fdnames;
                
                tensor1 = stm.eval();

                A    = tensor1.vals;
                tbl1 = tensor1.tbl;

                fds1  = tbl1.fdnames;
                ofds1 = fds1(~ismember(fds1, crossfds));
                
            else
                
                A    = prod.tensor1.vals;
                tbl1 = prod.tensor1.tbl;
                
            end

            if ~isempty(prod.reducefds2)
                
                assert(all(ismember(prod.reducefds2, ofds2)), ['There exist pre-reduced fields of tbl2 ' ...
                                                               'that do not belong to ofds2']);
                stm = SparseTensorMap(prod.tensor2);
                stm.toTbl    = projIndexArray(prod.tensor2.tbl, fds2(~ismember(fds2, prod.reducefds2)));
                stm.mergefds = stm.toTbl.fdnames;
                
                tensor2 = stm.eval();

                B    = tensor2.vals;
                tbl2 = tensor2.tbl;
                
                fds2  = tbl2.fdnames;
                ofds2 = fds2(~ismember(fds2, crossfds));
                
            else
                
                B    = prod.tensor2.vals;
                tbl2 = prod.tensor2.tbl;
                
            end
            

            lA = tbl1.gets(mergefds);
            iA = tbl1.gets(ofds1);
            jA = tbl1.gets(reducefds);

            lB = tbl2.gets(mergefds);
            jB = tbl2.gets(reducefds);
            kB = tbl2.gets(ofds2);

            if isempty(mergefds)
                mergefds = {''};
                lA = ones(size(lA, 1), 1, 'uint64');
                lB = ones(size(lB, 1), 1, 'uint64');
            end

            if isempty(reducefds)
                reducefds = {''};
                jA = ones(size(jA, 1), 1, 'uint64');
                jB = ones(size(jB, 1), 1, 'uint64');
            end

            if isempty(ofds1)
                ofds1 = {''};
                iA = ones(size(iA, 1), 1, 'uint64');
            end
            
            if isempty(ofds2)
                ofds2 = {''};
                kB = ones(size(kB, 1), 1, 'uint64');
            end

            
            % Fold the merge index l into the row and column index spaces of both
            % operands. Each nonzero of A then belongs to the column block l and each
            % nonzero of B to the row block l, so that the product sA*sB is block
            % diagonal in l : the cross-block entries are structurally zero and are
            % never computed. This gives the same multiplications as a loop over the
            % merge blocks, but with a single sparse product.

            [m_li, ~, i_A ] = unique([lA, iA], 'rows');
            [u_lj, ~, j_AB] = unique([[lA, jA]; [lB, jB]], 'rows');
            [m_lk, ~, k_B ] = unique([lB, kB], 'rows');

            j_A = j_AB(1 : size(jA, 1));
            j_B = j_AB(size(jA, 1) + (1 : size(jB, 1)));

            nI = size(m_li, 1);
            nJ = size(u_lj, 1);
            nK = size(m_lk, 1);

            if mrstVerbose > 0
                fprintf('sparse product of sizes %d x %d and %d x %d\n', nI, nJ, nJ, nK);
            end

            % Note that sparse() sums the values of the repeated indices, which is the
            % expected behaviour here.
            sA = sparse(i_A, j_A, A, nI, nJ);
            sB = sparse(j_B, k_B, B, nJ, nK);

            [i_C, k_C, C] = find(sA*sB);

            % The merge indices are given by the first columns of m_li (the ones of m_lk
            % coincide with them for each entry of the product, by construction).
            fdnames = {mergefds{:}, ofds1{:}, ofds2{:}};
            inds    = [m_li(i_C, :), m_lk(k_C, numel(mergefds) + 1 : end)];

            rmfds = cellfun(@(str) strcmp('', str), fdnames);

            fdnames = fdnames(~rmfds);
            inds    = inds(:, ~rmfds);
            
            tbl = IndexArray([]                , ...
                             'fdnames', fdnames, ...
                             'inds'   , inds);

            C = SparseTensor(C, tbl);

        end
        
    end
    
end

