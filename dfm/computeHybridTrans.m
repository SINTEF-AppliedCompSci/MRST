function [G,ft_hybrid] = computeHybridTrans(G,T)
% Computes the hybrid-hybrid transmissibilities between a hybrid cell 1 and 2
% using the star-delta transformation:
%
%  t_12 = t_1*t_2/sum(t_k)
%
%  for all cells k connected to an intersection.
%
% SYNOPSIS
%
%   [G, ft_hybrid] = computeHybridTrans(G,T)
%
% PARAMETERS:
%
%   G          - Grid data structure.
%   T          - Half-Transmissibilities one for each cell face.
%
% RETURNS:
%   G          - Grid structure with cell-to-cell connections added
%   ft_hybrid  - Transmissibilities for each hybrid cell2cell connection.
%
% Copyright 2011-2012 University of Bergen
%
% This file is licensed under the GNU General Public License v3.0.


% pick out the cellface index to the hybrid faces
[ind,loc] = ismember(G.cells.faces,G.hybridNeighbors.faces);

% pick out the hybrid half-Transmissibilities.
Th(loc(ind)) = T(ind);

% number of hybrid connections
nh = length(G.hybridNeighbors.facePos) - 1;

% mapping of the hybrid connections
hc = rldecode((1:nh)',diff(G.hybridNeighbors.facePos));

% compute the sum of the half transmissibilites related to each intersection.
sumTh = accumarray(hc, Th, [nh, 1]);

% compute the hybrid-hybrid transmissibility
ft_hybrid = prod(reshape(Th(G.hybridNeighbors.neighbors),sum(G.hybridNeighbors.n),2),2)./...
                                        rldecode(sumTh,G.hybridNeighbors.n);

% store cell2cell connections in the grid structure.
G.cells.neighbors = double(G.faces.neighbors(G.hybridNeighbors.faces(G.hybridNeighbors.neighbors)));
