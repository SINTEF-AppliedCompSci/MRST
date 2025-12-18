classdef SparseTensor
%
% SYNOPSIS:
%   SparseTensor(varargin)
%
% DESCRIPTION:
%
% Class for Sparse Tensor, which simply aggregate the value of the tensor with the corresponding IndexArray which gives the sparsity patten
%
% PARAMETERS:
%   - vals
%   - tbl
%
% RETURNS:
%   class instance
%
% EXAMPLE:
%   `computeVagTrans`
%
% SEE ALSO:
%   `IndexArray`, `TensorProd`.
    
    properties (SetAccess = immutable)
    
        vals % Value of the coefficients of the tensor
        tbl  % index array for the sparsity

    end
   
    methods
        
        function tensor = SparseTensor(vals, tbl)

            tensor.vals = vals;
            tensor.tbl  = tbl;

            assert(length(vals) == tbl.num, 'not matching values and index array');
            
        end

        function tensor = sparsify(tensor)

            % Remove zero values
            
            ind = (tensor.vals == 0);
            
            tensor.vals     = tensor.vals(~ind);
            tensor.tbl.inds = tensor.tbl.inds(~ind, :);
            
        end

        function tensor = setvals(tensors, vals, inds)

            upvals = tensor.values;
            upvals(inds) = vals;
            
            tensor = SparseTensor(upvals, A1.tbl);
    
        end
        
        function tensor = times(tensor, coef)

            tensor = SparseTensor((tensor.vals).*coef, tensor.tbl);

        end
        
        function v = values(tensor)
            
            v = tensor.vals;
            
        end

        function inds = indices(tensor)
            
            inds = tensor.tbl.inds;
            
        end
        
        function fds = fdnames(tensor)
            
            fds = tensor.tbl.fdnames;
            
        end
        
        function n = num(tensor)

            n = tensor.tbl.num;
            
        end

    end
    
end
