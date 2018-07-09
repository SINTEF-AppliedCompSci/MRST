function transmult = transmultEDFM(G,tol)
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

nummatcells=G.Matrix.cells.num;
numfraccells=G.cells.num-nummatcells;
matgeomnbors=sort(G.Matrix.faces.neighbors,2);
oldarealist=G.faces.areas;
newarealist=oldarealist;

% Get list of all matrix cells connected to each frac cell
for i=(1:numfraccells)+nummatcells
    % we search in the second column of G.nnc.cells because the fracture
    % cells are put in the second column in fracturematrixNNC3D.m
    index1= G.nnc.cells(:,2)==i;
    
    % pick only boundary fracture matrix connections. This will also help
    % filter out frac-frac connections
    index2=strcmp(G.nnc.type,'fracmat boundary');
    
    index=index1 & index2;
    
    % list of all matrix cells connected to frac cell i
    connectedmatcells=G.nnc.cells(index,1);
    
    matcellcombi=sort(nchoosek(connectedmatcells,2),2);
    
    lia=ismember(matgeomnbors,matcellcombi,'rows');
        
    if ~any(lia)
        continue; % skip to next iteration if lia is all zeros
    end
    
    index_globfaces=find(lia)'; % this tells us which areas to change in G.faces.area and G.faces.normal
    
    for j=index_globfaces        
        % extract fracture-matrix intersection area
        matcells=matgeomnbors(j,:);
        lia=ismember(G.nnc.cells,[matcells',[i;i]],'rows');
        fraccellarea=mean(G.nnc.area(lia));
        
        % extract current face area
        oldarea=newarealist(j);
        
        % calculate new face area
        newarea=oldarea-fraccellarea;
        if newarea<tol
            newarea=0;
        end
        
        % update face area for geometric neighbor
        newarealist(j)=newarea;
    end    
end

transmult=newarealist./oldarealist;



end
