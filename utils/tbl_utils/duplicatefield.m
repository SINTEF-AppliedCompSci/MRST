function newtbl = duplicatefield(tbl, fdcell)
%
%
% SYNOPSIS:
%   function newtbl = duplicatefield(tbl, fdcell)
%
% DESCRIPTION: Duplicate an index
%
% PARAMETERS:
%   tbl    - Index table
%   fdcell - Names of the index to duplicate with names of the duplicated
%   indices. The syntax is {'name', {'dupname1', 'dupname2'}}
%
% RETURNS:
%   newtbl - The resulting index table
%
% EXAMPLE:
%
% SEE ALSO: `setupTableMapping`
%

    
    oldfd = fdcell{1};
    fd1 = fdcell{2}{1};
    fd2 = fdcell{2}{2};

    [a, fds] = convertTableToArray(tbl, {oldfd});
    a = [a(:, 1), a];
    fds = fds(2 : end);
    fds = {fd1, fd2, fds{:}};
    
    newtbl = convertArrayToTable(a, fds);
    
end
