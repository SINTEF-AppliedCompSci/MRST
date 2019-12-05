function newtbl = replacefield(tbl, fieldpairs)
%
%
% SYNOPSIS:
%   function newtbl = replacefield(tbl, fieldpairs)
%
% DESCRIPTION: Utility function to replace field names in an index table
%
% PARAMETERS:
%   tbl        - Index table
%   fieldpairs - The syntax is either fieldpairs = {'A', 'B'} to replace name
%   'A' with 'B' or fieldpairs = {{'A1', 'B1'}, {'A2', 'B2'}} to replace
%   names 'A1' and 'A2' with 'B1' and 'B2', respectively
%
% RETURNS:
%   newtbl - Resulting index table
%
% EXAMPLE:
%
% SEE ALSO: `setupTableMapping`
%

    newtbl = tbl;
    if iscell(fieldpairs{1})
        for i = 1 : numel(fieldpairs)
            newtbl = replacethisfield(newtbl, fieldpairs{i});
        end
    else
        fieldpair = fieldpairs;
        newtbl = replacethisfield(newtbl, fieldpair);
    end
end

function tbl = replacethisfield(tbl, fieldpair)
    oldfield = fieldpair{1};
    newfield = fieldpair{2};
    tbl.(newfield) = tbl.(oldfield);
    tbl = rmfield(tbl, oldfield);
end