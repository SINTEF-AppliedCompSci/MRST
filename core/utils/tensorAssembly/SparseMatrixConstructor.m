classdef SparseMatrixConstructor

    properties (SetAccess = immutable)

        tensor 

    end

    properties

        fromTbl % IndexArray for source space (space of x in y = Mx)
        toTbl   % IndexArray for image space  (space of y in y = Mx)
        
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
        
        function smc = SparseTensorConstructor(tensor)
            
            smc = merge_options(smc, varargin{:}); 

            smc.tensor = tensor;
            
            if isempty(smc.reducefds)
                smc.reducefds = {};
            end
            
            if isempty(smc.reducefds1)
                smc.reducefds1 = {};
            end
            
            if isempty(smc.reducefds2)
                smc.reducefds2 = {};
            end

            if isempty(smc.mergefds)
                smc.mergefds = {};
            end            
            
        end

        function M = eval(smc)

            prod = TensorProd;
            prod.tbl1 = smc.tensor.tbl;
            prod.tbl2 = smc.fromTbl;
            prod.tbl3 = smc.toTbl;
            
            prod.reducefds   = smc.reducefds;
            prod.mergefds    = smc.mergefds;
            prod.replacefds1 = smc.replacefds1;
            prod.replacefds2 = smc.replacefds2;
            prod.replacefds3 = smc.replacefds3;
            prod.reducefds1  = smc.reducefds1;
            prod.reducefds2  = smc.reducefds2;

            M = prod.setupMatrix(smc.tensor.values);
            
        end
        
    end
    
end

