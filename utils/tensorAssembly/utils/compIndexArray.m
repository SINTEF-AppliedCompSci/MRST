function iseq = compIndexArray(tbl1, tbl2)
%
%
% SYNOPSIS:
%   function iseq = compIndexArray(tbl1, tbl2)
%
% DESCRIPTION: Compares the two IndexArrays tbl1 and tbl2 (field names and indices)
%
% PARAMETERS:
%   tbl1 - IndexArray
%   tbl2 - IndexArray
%
% RETURNS:
%   iseq - true if equal false otherwise
%
% EXAMPLE:
%
% SEE ALSO:
%

    
    fds1 = tbl1.fdnames;
    fds2 = tbl2.fdnames;    
    
    if tbl1.isvirtual | tbl2.isvirtual
        % the comparison cannot be done since the Index arrays are virtual
        iseq = NaN;
        return
    end
  
    iseq = (numel(fds1) == numel(fds2));
    
    if iseq
        
        [~, fdind2] = ismember(fds1, fds2);
        
        mat1 = sortrows(tbl1.inds);
        mat2 = sortrows(tbl1.inds(:, fdind2));
        
        iseq = all(all((mat1 - mat2) == 0));
    end
    
end
