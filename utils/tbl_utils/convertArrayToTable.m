function tbl = convertArrayToTable(A, fds)
%
%
% SYNOPSIS:
%   function tbl = convertArrayToTable(A, fds)
%
% DESCRIPTION: Convert an array to an indexing table. The indexing table are
% given as structure with field names. Each field contains a vector of
% indices. All these vectors have same length. 
%
% PARAMETERS:
%   A   - Array
%   fds - Field names given to each of the column of A
%
% RETURNS:
%   tbl - The corresponding indexing table
%
% EXAMPLE:
%
% SEE ALSO: `convertTableToArray`, `setupTableMapping`
%

    sz = size(A);
    nfds = numel(fds);
    
    for ifield = 1 : nfds
        tbl.(fds{ifield}) = A(:, ifield);
    end
    tbl.num = sz(1);
end
