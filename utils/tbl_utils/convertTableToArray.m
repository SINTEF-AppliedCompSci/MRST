function [A, fds] = convertTableToArray(tbl, fdsfirst)
%
%
% SYNOPSIS:
%   function [A, fds] = convertTableToArray(tbl, fdsfirst)
%
% DESCRIPTION: Convert an indexing table to an array. The indexing table are
% given as structure with field names. Each field contains a vector of
% indices. All these vectors have same length.
%
% PARAMETERS:
%   tbl      - Indexing table
%   fdsfirst - We give the fields that we want to appear first and in the
%   given order in the array
%
% RETURNS:
%   A   - Corresponding array
%   fds - returns the name of the fields corresponding to the columns in A
%
% EXAMPLE:
%
% SEE ALSO: `convertArrayToTable`, `setupTableMapping`
%
    
    fds = fieldnames(tbl);
    
    fdsToRemove = {fdsfirst{:}, 'ind', 'num'};
    for ifield = 1 : numel(fdsToRemove)
        fds = fds(~strcmp(fds, fdsToRemove{ifield}));
    end
    
    fds = {fdsfirst{:}, fds{:}};
    A = [];
    A = tbl.(fds{1});
    for ifield = 2 : numel(fds)
        A = [A, tbl.(fds{ifield})];
    end
     
end
