function projtbl = projTable(tbl, fds)
%
%
% SYNOPSIS:
%   function projtbl = projTable(tbl, fds)
%
% DESCRIPTION: Project an index table on a subset of indices, meaning that
% given a table tbl with field tbl.A and tbl.B. If we project on the index
% name 'A', we construct a projtbl with field projtbl.A such that, for all i,
% there exists a unique j such that tbl.A(i) = projtbl.A(j) 
%
% PARAMETERS:
%   tbl - Index table
%   fds - field names along which we project. Syntax is {'A1', 'A2'} for example
%
% RETURNS:
%   projtbl - resulting index table
%
% EXAMPLE:
%
% SEE ALSO: `setupTableMapping`
%

    a = convertTableToArray(tbl, fds);
    nfds = numel(fds);
    a = a(:, (1 : nfds));
    a = unique(a, 'rows');
    projtbl = convertArrayToTable(a, fds);
end