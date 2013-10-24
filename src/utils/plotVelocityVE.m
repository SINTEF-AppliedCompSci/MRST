function plotVelocityVE(g, state)
%Plot two-dimensional velocity vectors in solution.
%
% SYNOPSIS:
%   plotVelocityVE(G, state)
%
% PARAMETERS:
%   G     - A grid structure.
%
%   state - An MRST solution structure containing a valid flux field.
%
% RETURNS:
%   Nothing, but displays a QUIVER plot of the lateral velocity components.
%
% SEE ALSO:
%   quiver, grid_structure.

%{
#COPYRIGHT#
%}

cellflux     = faceFlux2cellFlux(g, state.flux);

cellvelocity = cellFlux2cellVelocity(g, cellflux);

cellvelocity = bsxfun(@rdivide, cellvelocity, ...
                      sqrt(sum(cellvelocity .^ 2, 2)));

quiver(g.cells.centroids(:,1), ...
       g.cells.centroids(:,2), ...
       cellvelocity(:,1), cellvelocity(:,2));
end

function cellvel = cellFlux2cellVelocity(G, cellflux)

   hf=G.cells.faces(:,1);
   cells = rldecode([1:G.cells.num]',diff(G.cells.facePos));
   hfc = G.faces.centroids(hf,:)-G.cells.centroids(cells,:);
   tmp = bsxfun(@times,hfc,cellflux)';
   cellvel(:,1)=accumarray(cells,tmp(1,:));
   cellvel(:,2)=accumarray(cells,tmp(2,:));
   cellvel = bsxfun(@rdivide,cellvel,G.cells.volumes);

end
