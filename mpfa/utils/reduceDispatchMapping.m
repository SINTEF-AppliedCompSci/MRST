function map = reduceDispatchMapping(tbl1, tbl2, crossfield, varargin)

    ind1 = tbl1.(crossfield);
    ind2 = tbl2.(crossfield);
    n1 = tbl1.num; 
    n2 = tbl2.num;
    n = max(n1, n2);
    
    imap1 = sparse(ind1, (1 : n1)', 1, n, n1);
    imap2 = sparse(ind2, (1 : n2)', 1, n, n2);    
    map = imap2'*imap1;
end
