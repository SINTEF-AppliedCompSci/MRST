classdef SparseTensorGroupBy < handle

    properties (SetAccess = immutable)

        vals
        
    end

    properties (SetAccess = private)

        iagb % current index in the group list
        prod % prod in ofdnames to setup matrix
        
    end

    properties
        
        replacefds1
        fds2
        replacefds2
        fds3
        replacefds3
        
    end
    
    properties (Dependent)

        i
        
    end
    
    methods

        function stgb = SparseTensorGroupBy(tensor, gfdnames)

            stgb.iagb = IndexArrayGroupBy(tensor.tbl, gfdnames);
            stgb.vals = tensor.vals;
            stgb.i = 0;

        end

        function gt = gtensor(stgb)

            gt = SparseTensor(stgb.vals(stgb.iagb.startInd() : stgb.iagb.endInd()), stgb.iagb.grouptbl());
            
        end

        function setupprod(stgb)
            
            gt = stgb.gtensor();

            prod = TensorProd();
            
            prod.tbl1 = gt.tbl;
            prod.replacefds1 = stgb.replacefds1;
            
            prod.tbl2 = gt.tbl.proj(stgb.fds2);
            prod.replacefds2 = stgb.replacefds2;
                
            prod.tbl3 = gt.tbl.proj(stgb.fds3);
            prod.replacefds3 = stgb.replacefds3;

            prod.reducefds = stgb.fds2;
            
            stgb.prod = prod.setup();

        end
        
        function M = gmatrix(stgb)

            gt = stgb.gtensor();
            
            stgb.setupprod();
            
            M = stgb.prod.setupMatrix(gt.values);
            
        end
        
        function reset(stgb)
            
            stgb.i = 0;
            
        end
        
        function set.i(stgb, val)
            
            stgb.iagb.i = val;
            
        end

        function val = get.i(stgb)
            
            val = stgb.iagb.i;
            
        end
        
    end

    methods(Static)

        function fds = replacefds(fds, replacements)
            for irep = 1 : numel(replacements)
                fd1 = replacements{irep}{1};
                fd2 = replacements{irep}{2};
                [isok, lia] = ismember(fd1, fds);
                assert(isok, 'fieldname not found');
                fds{lia} = fd2;
            end
        end
        
    end
    
end

