function tbl = addLocInd(tbl, locindname)
%
%
% SYNOPSIS:
%   function tbl = addLocInd(tbl, locindname)
%
% DESCRIPTION: Add a new field with a local indexing (1,2, ..., size of
% table) to the given indexing table.
%
% PARAMETERS:
%   tbl        - Indexing table
%   locindname - Chosen name for the local index
%
% RETURNS:
%   tbl - Indexing table with new field giving a local index (from 1 to size
%   of table)
%
% EXAMPLE:
%
% SEE ALSO:
%

    tbl.(locindname) = (1 : tbl.num)';
end
