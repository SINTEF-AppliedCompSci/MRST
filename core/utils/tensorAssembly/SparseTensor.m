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
    
    properties
    
        vals % Value of the coefficients of the tensor
        tbl  % index array for the sparsity

    end
   
    methods
        
        function tensor = SparseTensor(vals, tbl)

            tensor.vals = vals;
            tensor.tbl  = tbl;
            
        end

    end
    
end
