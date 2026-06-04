function [iseq, tens1, tens2] = compareSparseTensor(tens1, tens2)

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

    if all(tens1.vals == tens2.vals)
        iseq = true;
    else
        iseq = false;
    end
    
    
end
