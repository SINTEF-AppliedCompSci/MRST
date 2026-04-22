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

    
    properties
    
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

        
        function tensor = times(tensor, coef)

            tensor.vals = (tensor.vals).*coef;

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

        function tensor2 = sort(tensor, fdnames)

            [tbl, dispind] = sortIndexArray(tensor.tbl, fdnames, 'keepAllFields', true);

            tensor2 = tensor;
            
            tensor2.vals = tensor.vals(dispind);
            tensor2.tbl  = tbl;
            
        end
        
        function table = print(tensor, varargin)

            tbl = tensor.tbl;
            
            opt = struct('range', (1 : tbl.num)', ...
                         'fdnames', {tbl.fdnames});

            opt = merge_options(opt, varargin{:});
            
            [found, indcol] = ismember(opt.fdnames, tbl.fdnames);
            assert(all(found), 'some field names were not found');

            indrow = opt.range;
            inds = tbl.inds(indrow, indcol);

            vals = [tensor.vals(indrow), double(inds)];
            

            fdnames = horzcat({'values'}, tbl.fdnames(indcol));
            fmt = sprintf('%%10s %s \n', strcat(repmat('%10s ', 1, numel(indcol))));
            fprintf(fmt, fdnames{:});
            fprintf('\n');
            fmt = sprintf('%%10g %s \n', strcat(repmat('%10d ', 1, numel(indcol))));
            fprintf(fmt, vals');
            
        end

    end
    
end
