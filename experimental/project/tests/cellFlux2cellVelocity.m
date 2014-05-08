function cellvel = cellFlux2cellVelocity(G, cellflux)
%Undocumented internal helper in Vertical Equilibrium module.

%{
#COPYRIGHT#
%}

   %ff = mcolon(G.cells.facePos(1:end-1),G.cells.facePos(2:end));
   %hf = G.cells.faces(G.cells.facePos(1:end-1),G.cells.facePos(2:end),1);
   hf=G.cells.faces(:,1);
   cells = rldecode([1:G.cells.num]',diff(G.cells.facePos));
   hfc = G.faces.centroids(hf,:)-G.cells.centroids(cells,:);
   tmp = bsxfun(@times,hfc,cellflux)';
   %cellvelpos = [repmat(cells,size(tmp,1),1),repmat([1:size(tmp,1)]',size(tmp,2),1)];
   %cellvel=accumarray(cellvelpos,tmp(:));
   cellvel(:,1)=accumarray(cells,tmp(1,:));
   cellvel(:,2)=accumarray(cells,tmp(2,:));
   cellvel = bsxfun(@rdivide,cellvel,G.cells.volumes);

end
