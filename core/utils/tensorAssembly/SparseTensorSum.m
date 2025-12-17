classdef SparseTensorSum

    properties (SetAccess = immutable)

        tensor1 % first tensor argument
        tensor2 % second tensor argument

    end

    properties

        replacefds1  % Possibility to change field names of first table
                     % (before setting up product)
        replacefds2  % Possibility to change field names of second table
                     % (before setting up product)

        coefs1
        coefs2
        
    end
   
    methods
        
        function tsum = SparseTensorSum(tensor1, tensor2, varargin)
            
            opts = struct('replacefds1', [], ...
                          'replacefds2', [], ...
                          'coefs1'     , [], ...
                          'coefs2'     , []);
            
            tsum = merge_options(tsum, varargin{:}); 

            tsum.tensor1 = tensor1;
            tsum.tensor2 = tensor2;
            
        end

        function C = eval(tsum)

            tblA = tsum.tensor1.tbl;
            tblB = tsum.tensor2.tbl;

            if ~isempty(tsum.replacefds1)
                tblA = replacefield(tblA, tsum.replacefds1);
            end
            
            if ~isempty(tsum.replacefds2)
                tblB = replacefield(tblB, tsum.replacefds2);
            end

            fdnames = tblA.fdnames;
            fdsB    = tblB.fdnames;

            assert(numel(fdsB) == numel(fdnames), 'not matching field names for the two Tensors');

            [ok, inds] = ismember(fdsB, fdnames);

            assert(all(ok), 'not matching field names for the two Tensors');

            iA = tblA.inds;
            iB = tblB.inds;
            iB = iB(:, inds);
            
            [m_b_i, ~, b_i] = unique([iA; iB], 'rows');
            b_iA = b_i(1 : size(iA, 1));
            b_iB = b_i(size(iA, 1) + (1 : size(iB, 1)));

            A = sparse(b_iA, ones(size(b_iA)), tsum.tensor1.values, size(m_b_i, 1), 1);
            B = sparse(b_iB, ones(size(b_iB)), tsum.tensor2.values, size(m_b_i, 1), 1);

            if ~isempty(tsum.coefs1)
                A = tsum.coefs1.*A;
            end
            
            if ~isempty(tsum.coefs2)
                B = tsum.coefs2.*B;
            end
            
            [i, j, v] = find(A + B);
            
            tbl = IndexArray([]                , ...
                             'fdnames', fdnames, ...
                             'inds'   , m_b_i(i, :));

            C = SparseTensor(v, tbl);
            
        end
        
    end
    
end

