function tbl = sortTable(tbl, fds, varargin);
%
%
% SYNOPSIS:
%   function tbl = sortTable(tbl, fds, varargin)
%
% DESCRIPTION: Sort in increasing order the rows of indices given by the
% indexing table tbl, following the order given by fds. (We use sortrows)
%
% PARAMETERS:
%   tbl      - Indexing table
%   fds      - Field names which give the order of the columns for the row
%              indexing sort
%   varargin - 
%
% KEYWORD ARGUMENTS:
%
%   keepAllFields - We keep all the fields for the table, even those that are
%                   not given in fds but belong to tbl. The default value is false. This
%                   option has not been properly tested and we recommend not using it for the moment.
%   
% RETURNS:
%   tbl - Indexing table that has been sorted by row.
%
% EXAMPLE:
%
% SEE ALSO: `setupTableMapping`
%

    opt = struct('keepAllFields', false);
    opt = merge_options(opt, varargin{:});
    
    a = convertTableToArray(tbl, fds);
    a = sortrows(a);
    
    if opt.keepAllFields
        warning('option not really tested yet!');
        ofds = fieldnames(tbl);
    
        fdsToRemove = {fds{:}, 'num'};
        for ifield = 1 : numel(fdsToRemove)
            ofds = ofds(~strcmp(ofds, fdsToRemove{ifield}));
        end
        
        fds = {fds{:}, ofds{:}};
    end
    
    tbl = convertArrayToTable(a, fds);
end
