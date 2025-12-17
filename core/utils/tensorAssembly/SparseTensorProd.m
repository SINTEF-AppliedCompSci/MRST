classdef SparseTensorProd

    properties (SetAccess = immutable)

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
            
            opts = struct('replacefds1', [], ...
                          'replacefds2', [], ...
                          'replacefds3', [], ...
                          'reducefds'  , [], ...
                          'reducefds1' , [], ...
                          'reducefds2' , [], ...
                          'mergefds'   , []);
            
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
            
            tbl1 = prod.tensor1.tbl;
            tbl2 = prod.tensor2.tbl;

            if ~isempty(prod.replacefds1)
                tbl1 = replacefield(tbl1, prod.replacefds1);
            end
            if ~isempty(prod.replacefds2)
                tbl2 = replacefield(tbl2, prod.replacefds2);
            end

            % some sanity checks on the table's field names.
            fds1 = tbl1.fdnames;
            fds2 = tbl2.fdnames;

            assert(all(ismember(mergefds, fds1)), ['There exist merge fields that do ' ...
                                'not belong to fields of first table']);
            assert(all(ismember(mergefds, fds2)), ['There exist merge fields that do ' ...
                                'not belong to fields of second table']);
            ofds1 = fds1(~ismember(fds1, crossfds));
            ofds2 = fds2(~ismember(fds2, crossfds));
            assert(all(~ismember(ofds1, ofds2)) & all(~ismember(ofds2, ofds1)), ...
                   ['There exist fields with same name in first and second ' ...
                    'table that are neither merged or reduced.']);
            
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

            A = prod.tensor1.vals;
            B = prod.tensor2.vals;
            
            [m_l, ~, b_l] = unique([lA; lB], 'rows');
            b_lA = b_l(1 : size(lA, 1));
            b_lB = b_l(size(lA, 1) + (1 : size(lB, 1)));

            [~, ~, b_j] = unique([jA; jB], 'rows');
            b_jA = b_j(1 : size(jA, 1));
            b_jB = b_j(size(jA, 1) + (1 : size(jB, 1)));

            [m_i, ~, b_iA] = unique(iA, 'rows');
            [m_k, ~, b_kB] = unique(kB, 'rows');

            b_l = intersect(b_lA, b_lB, 'sorted');

            LIA = ismember(b_lA, b_l);
            b_lA = b_lA(LIA);
            b_iA = b_iA(LIA);
            b_jA = b_jA(LIA);
            A    = A(LIA);

            LIA = ismember(b_lB, b_l);
            b_lB = b_lB(LIA);
            b_jB = b_jB(LIA);
            b_kB = b_kB(LIA);
            B    = B(LIA);

            [b_lA, isort] = sort(b_lA);
            b_iA = b_iA(isort);
            b_jA = b_jA(isort);
            A = A(isort);
 
            [b_lB, isort] = sort(b_lB);
            b_jB = b_jB(isort);
            b_kB = b_kB(isort);
            B = B(isort);

            b_lA_b_iA = unique([b_lA, b_iA], 'rows');
            [~, nI] = rlencode(b_lA_b_iA(: , 1));
            
            b_lB_b_kB = unique([b_lB, b_kB], 'rows');
            [~, nK] = rlencode(b_lB_b_kB(: , 1));

            nalloc = sum(nI.*nK);

            [u_b_lA, n_lA] = rlencode(b_lA);
            [u_b_lB, n_lB] = rlencode(b_lB);

            assert(all(u_b_lA == u_b_lB));
            u_b_l = u_b_lA;
            
            startA   = 1;
            endA     = n_lA(1);
            
            startB   = 1;
            endB     = n_lB(1);
            
            b_iC_all = zeros(nalloc, 1, 'uint64');
            b_kC_all = zeros(nalloc, 1, 'uint64');
            C_all    = zeros(nalloc, 1, 'double');
            
            n_u_b_l = [];
            
            posC = 0;

            if mrstVerbose > 0
                fprintf('number of sparse multiplication: %d\n ', numel(u_b_l));
            end

            for ind_u_b_l = 1 : numel(u_b_l)
                
                ind = (startA : endA);
                sA = sparse(b_iA(ind), b_jA(ind), A(ind));

                ind = (startB : endB);
                sB = sparse(b_jB(ind), b_kB(ind), B(ind));

                sC = sA * sB;

                [b_iC, b_kC, C] = find(sC);                    

                nC = size(C, 1);
                
                n_u_b_l(end + 1) = nC;
                
                b_iC_all(posC  + (1 : nC)) =  b_iC;
                b_kC_all(posC  + (1 : nC)) =  b_kC;
                C_all(posC  + (1 : nC))    =  C ;
                
                posC = posC + nC;

                if ind_u_b_l < numel(u_b_l)
                    startA = endA + 1;
                    endA   = startA + n_lA(ind_u_b_l + 1) - 1;                        
                    startB = endB + 1;
                    endB   = startB + n_lB(ind_u_b_l + 1) - 1;
                end

            end

            if nalloc < posC
                warning('allocation was not enough');
                fprintf('nalloc : %d, posC : %d, difference : %d\n', nalloc, posC, abs(nalloc - posC));
            end
            
            C = C_all(1 : posC);

            lC_all = rldecode(m_l(u_b_l, :), n_u_b_l);
            
            fdnames = {mergefds{:}, ofds1{:}, ofds2{:}};
            inds    =  [lC_all, m_i(b_iC_all(1 : posC), :), m_k(b_kC_all(1 : posC), :)];

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

