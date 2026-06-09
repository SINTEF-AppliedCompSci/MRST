function [iseq, tens1, tens2] = compareSparseTensor(tens1, tens2, options)

    if nargin < 3
        options = [];
    end

    options = setDefaultStructField(options, 'useTolerance', false);
    options = setDefaultStructField(options, 'toleranceType', 'relative');
    options = setDefaultStructField(options, 'tolerance', 1e-12);
    
    fdnames1 = tens1.tbl.fdnames;
    fdnames2 = tens2.tbl.fdnames;
    
    if ~(all(ismember(fdnames1, fdnames2)) && all(ismember(fdnames1, fdnames2)))
        iseq = false
        return
    end
    
    tbl = concatIndexArray(tens1.tbl, tens2.tbl, {}, 'checkUnique', false, 'removeDuplicates', true);

    stm = SparseTensorMap(tens1);
    stm.toTbl = tbl;
    stm.mergefds = 'allTo';

    tens1 = stm.eval();

    stm = SparseTensorMap(tens2);
    stm.toTbl = tbl;
    stm.mergefds = 'allTo';

    tens2 = stm.eval();

    if options.useTolerance

        tol = options.tolerance;
        
        switch options.toleranceType

          case 'relative'
            m = max(abs(tens1.vals), abs(tens2.vals));
            if all(abs(tens1.vals(m > 0) - tens2.vals(m > 0))./m(m > 0) < tol)
                iseq = true;
            else
                iseq = false;
            end
            
          case 'absolute'
            
            if all(abs(tens1.vals - tens2.vals) < tol)
                iseq = true;
            else
                iseq = false;
            end
            
          otherwise
            
            error('toleranceType not recognized')
            
        end
        
    else
        
        if all(tens1.vals == tens2.vals)
            iseq = true;
        else
            iseq = false;
        end
    end
    
end
