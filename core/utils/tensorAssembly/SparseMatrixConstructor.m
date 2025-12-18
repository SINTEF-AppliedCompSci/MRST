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
        
        function gen = SparseTensorConstructor(tensor)
            
            opts = struct('replacefds1', [], ...
                          'replacefds2', [], ...
                          'replacefds3', [], ...
                          'reducefds'  , [], ...
                          'reducefds1' , [], ...
                          'reducefds2' , [], ...
                          'mergefds'   , []);
            
            gen = merge_options(gen, varargin{:}); 

            gen.tensor = tensor;
            
            if isempty(gen.reducefds)
                gen.reducefds = {};
            end
            
            if isempty(gen.reducefds1)
                gen.reducefds1 = {};
            end
            
            if isempty(gen.reducefds2)
                gen.reducefds2 = {};
            end

            if isempty(gen.mergefds)
                gen.mergefds = {};
            end            
            
        end

        function M = eval(gen)

            prod = TensorProd;
            prod.tbl1 = gen.tensor.tbl;
            prod.tbl2 = gen.fromTbl;
            prod.tbl3 = gen.toTbl;
            
            prod.reducefds   = gen.reducefds;
            prod.mergefds    = gen.mergefds;
            prod.replacefds1 = gen.replacefds1;
            prod.replacefds2 = gen.replacefds2;
            prod.replacefds3 = gen.replacefds3;
            prod.reducefds1  = gen.reducefds1;
            prod.reducefds2  = gen.reducefds2;

            M = prod.setupMatrix(gen.tensor.values);
            
        end
        
    end
    
end

