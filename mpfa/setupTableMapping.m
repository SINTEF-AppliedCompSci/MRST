function [op, colind, rowind] = setupTableMapping(var1, var2, varargin)
    
    fds = varargin;
    nfds = numel(fds);
    
    for ifield = 1 : nfds
        fieldname = fds{ifield};
        
        ind1 = var1.(fieldname);
        ind2 = var2.(fieldname);
        
        n1 = numel(ind1);
        n2 = numel(ind2);
        n = max([ind1; ind2]);
        
        iop1 = sparse(ind1, (1 : n1)', 1, n, n1);
        iop2 = sparse(ind2, (1 : n2)', 1, n, n2);
        
        iop = iop2'*iop1;
        if ifield == 1
            op = iop;
        else
            op = op.*iop;
        end
    end
    
    [colind, rowind] = find(op);
    
end
