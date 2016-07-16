function [ncoords, ixs] = uniqueNodes(coords)
% Extract from the provided list of (non-unique) coordinates a unique set,
% and provide the corresponding indexing.
%
% SYNOPSIS:
%   function [ncoords, ixs] = uniqueNodes(coords)
%
% DESCRIPTION:
%
% PARAMETERS:
%   coords - N x D matrix of the D-dimensional coordinates for a set of
%            points.  Each row represents one point.
%
% RETURNS:
%   ncoords - M x D matrix of M points with unique nodes, resulting from
%             removing all duplicates from the 'coords' matrix.
%   ixs     - Indexing.  The point corresponding to row i in the original
%             'coords' matrix will find its coordinates at row ixs(i) in the
%             returned matrix 'ncoords'.
%

   [coords, ix] = sortrows(coords);
   
   % deciding which coordinates to keep (discarding duplicates)
   keep = [true; any(diff(coords), 2)];
   
   % Determining number of duplicates for each coord
   [~, dup] = rlencode(keep);
   
   ncoords = coords(keep,:);
   
   % Determining indexing into vector without duplicates
   ixs = zeros(size(coords, 1), 1);
   ixs(keep) = [1:sum(keep)]';
   
   diff_keep = diff(keep);
   removed_start_ixs = find(diff_keep == -1) + 1;
   removed_end_ixs   = find(diff_keep ==  1);
   if numel(removed_end_ixs) < numel(removed_start_ixs)
      removed_end_ixs = [removed_end_ixs; numel(keep)];
      assert(numel(removed_end_ixs) == numel(removed_start_ixs));
   end
   
   lengths = removed_end_ixs - removed_start_ixs + 1;
   vals = ixs(removed_start_ixs-1);
   
   ixs(~keep) = rldecode(vals, lengths);

   [~, ix] = sort(ix);
   ixs = ixs(ix);
   
end
