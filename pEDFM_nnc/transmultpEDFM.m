function transmultp = transmultpEDFM(G,tol)
% This function calculates transmissibility multipliers for geometric
% neighbors. For fracture cells that reside on the boundary of two
% neighboring cells, since a transmissibility was established between the
% frac cell and the matrix cells, the transmissibility between the matrix
% cells need to be reduced according to the area of the intersecting
% fracture cell.
%
% The multipliers can be multiplied with the full transmissibility of
% neighboring cells to correct this behavior.

% METHODOLOGY: Go through every global fracture cell and gather all the
% matrix cells that they are connected to. Filter out the connections that
% are not boundary types. Once we have this list, then for any two matrix cells
% connected to a fracture cell, if they are also geometric neighbours, then
% reduce the area according to data from G.nnc.area. Do this for all
% combinations of matrix cells for all fracture cells. Use the new area
% calculated and the original face area to generate multipliers.

% nummatcells=G.Matrix.cells.num;
oldarealist=G.faces.areas;
newarealist=oldarealist;


numMp_M = numel(G.nnc.pMMneighs(:,2));
pMFidxInNNC = ( (1:numMp_M) + numel(G.nnc.normal(:,1)) )';
projectedArea = G.nnc.area(pMFidxInNNC);
globalFaceIdx = zeros(numMp_M,1);
for j=1:numMp_M
    lia2 = ismember(G.Matrix.faces.neighbors,G.nnc.pMMneighs(j,:),'rows');
    globalFaceIdx(j) = find(lia2);
end  
oldArea = newarealist(globalFaceIdx);
newArea = oldArea - projectedArea;
newArea(newArea<tol) = 0;
newarealist(globalFaceIdx) = newArea;

transmultp = newarealist./oldarealist;



%%This approach works when there is no fracture intersection
% [pMMneighsUniq, ia, ic] = unique(G.nnc.pMMneighs,'rows','stable');  %pMMneighsUniq = G.nnc.pMMneighs(ia,:);
% numUniqMp_M = numel(pMMneighsUniq(:,2));
% % lia = ismember(G.Matrix.faces.neighbors,pMMneighsUniq,'rows');
% % globalFaceIdx = find(lia); %this sorts the indices and messes up the result
% pMFidxInNNC = ( (1:numel(G.nnc.pMMneighs(:,2))) + numel(G.nnc.normal(:,1)) )';
% areaAll = G.nnc.area(pMFidxInNNC); 
% projectedArea = areaAll(ia); %these are the areas for the unique pairs of pM-M cells 
% 
% globalFaceIdx = zeros(numUniqMp_M,1);
% for j=1:numUniqMp_M
%     lia2 = ismember(G.Matrix.faces.neighbors,pMMneighsUniq(j,:),'rows');
%     globalFaceIdx(j) = find(lia2);
% end    
% 
% 
% oldArea = newarealist(globalFaceIdx);
% newArea = oldArea - projectedArea;
% 
% % test = newArea(newArea<tol);
% % numel(newArea(newArea<tol))
% 
% newArea(newArea<tol) = 0;
% newarealist(globalFaceIdx) = newArea;
% 
% transmultp = newarealist./oldarealist;



%% This approach fails because find(lia) messes up with the ordering of the globalFaceIdx
% lia = ismember(G.Matrix.faces.neighbors,G.nnc.pMMneighs,'rows');
% globalFaceIdx = find(lia);
% %Recall that G.nnc.pMMneighs only lists pMM cell pairs, and is a smaller 
% %perfect subset of G.nnc.cells. Also, G.nnc.normal only lists the
% %frac-related cell pairs, and is a perfect subset of G.nnc.cells.
% %Loosely speaking, in terms of the neighbors being referred to 
% %(and not the actual content of these variables, we can say that
% % (1) G.nnc.cells = G.nnc.normal "UNION" G.nnc.pMMneighs
% % (2) G.nnc.normal "INTERSECTION" G.nnc.pMMneighs = "NULL SET"
% 
% %The line below then gives us the rows of G.nnc.cells that corresponds to
% %the rows in the G.nnc.pMMneighs
% %Last "numel(G.nnc.normal(:,1)) + 1:numel(globalFaceIdx)" in G.nnc.cells
% %are the indices of the pM-F ("fracpmat")connections
% pMFidxInNNC = ( (1:numel(globalFaceIdx)) + numel(G.nnc.normal(:,1)) )';
% projectedArea = G.nnc.area(pMFidxInNNC);
% oldArea = newarealist(globalFaceIdx);
% newArea = oldArea - projectedArea;
% newArea(newArea<tol) = 0;
% newarealist(globalFaceIdx) = newArea;
% 
% transmultp = newarealist./oldarealist;



end