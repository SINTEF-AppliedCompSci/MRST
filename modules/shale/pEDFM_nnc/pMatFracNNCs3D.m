function G = pMatFracNNCs3D(G,tol, varargin)
    opt = struct('defaultFace', 'bottom');
    opt = merge_options(opt, varargin{:});
    
%     startIdx = numel(G.nnc.cells(:,1))+1; %this will be the index of first projected cell
    idx = strcmp(G.nnc.type,'fracmat interior');
    mcells = G.nnc.cells(idx,:); %Nx2 array of matrix and fracture cell pairs (that are within a matrix cell)

    faceStartIdx = G.cells.facePos(mcells(:,1));   %indices of the first face of each matrix cell in mcells(:,1)
    faceEndIdx = G.cells.facePos(mcells(:,1)+1)-1; %indices of the last face of each matrix cell in mcells(:,1)

    %assuming 3d
    numInteriorFcells = numel(mcells(:,2));
    nRows4pNNC = numInteriorFcells*3;
    p.MFcells = zeros(nRows4pNNC,2);
    p.MMneighs = zeros(nRows4pNNC,2);
    p.Type = cell(nRows4pNNC,1);
    [p.Type{1:nRows4pNNC}] = deal('ProjMatrix outside');
    p.Trans = zeros(nRows4pNNC,1);
    p.Area = zeros(nRows4pNNC,1);
    p.Dir = zeros(nRows4pNNC,1);
    
    fracMatAreas = G.nnc.area(idx,:);
    tic
    %Loop through each matrix cell that hosts a fracture and find its
    %projected matrix cell neighbors (called pmcells)
    for i=1:numInteriorFcells
        faceIdxs_i = G.cells.faces(faceStartIdx(i):faceEndIdx(i)); %indices of faces that bound mcells(i,1)
        faceSelected = zeros(numel(faceIdxs_i),1);

        projectionCells = G.faces.neighbors(faceIdxs_i,:); %list of mcells(i,1) and it's neighbors

        faceCentroids = G.faces.centroids(faceIdxs_i,:);
        %TO DO: remember to check for boundary cells that host a fracture.
        %Neighbor matrix cell will be 0 in this case, and
        %numel(projectionCells(:,1) may or may not be less than 6

        %In projectionCells, row1 is Left,Me; row2 is Me,Right; row3 is
        %Front,Me; row4 is Me,Back; row5 is Bottom,Me; row6 is Me,Top
        fraccellcentroid = G.cells.centroids(mcells(i,2),:); %centroid of current frac cell


        %compute distance between frac centroid and face centroid in all 6 directions 
        distLeft = norm(faceCentroids(1,:)-fraccellcentroid);
        distRight = norm(faceCentroids(2,:)-fraccellcentroid);
        distFront = norm(faceCentroids(3,:)-fraccellcentroid);
        distBack = norm(faceCentroids(4,:)-fraccellcentroid);
        distBottom = norm(faceCentroids(5,:)-fraccellcentroid);
        distTop = norm(faceCentroids(6,:)-fraccellcentroid);
        
        isEquiDist = zeros(3,1); %initialize isEquiDist to 0 in all 3 directions
        projectionCell = zeros(3,1); %initialize selected projection cell in all 3 directions

        %compute distance between frac centroid and pmcell centroid in x, y & z directions (dist_MpF_x,y,z)     
        if (distLeft+tol)<distRight    %distLeft<distRight
            %left face is selected
%             faceSelected(1) = 1; 
%             localIdx = 3*(i-1)+1;
%             projectionCell(1) = projectionCells(1,1);
%             %projection cells always in first column and frac cell in second
%             %column because matrix cells come before frac cells in G.cells
%             pMFcells(localIdx,:) = [projectionCell(1), mcells(i,2)];
%             %algorithm for pneighs is such that you enumerate projected matrix
%             %in first column of pneighs if left, front and bottom cells are 
%             %the selected projection cells. The matrix which hosts the frac cell
%             %is put in the second column in these cases. This is reversed if 
%             %the projection cells aree the right, back and top cells.
%             pMMneighs(localIdx,:) = [projectionCell(1), mcells(i,1)];
%             [pType{localIdx}, pArea(localIdx), pTrans(localIdx)] = ...
%                     computeNNCprojArea_n_Trans(G,pMMneighs(localIdx,:),pMFcells(localIdx,:),mcells(i,:),fracMatAreas(i));
            %left
            currFaceID = 1;
            faceSelected(currFaceID) = 1;                        
            localIdx = 3*(i-1)+ceil(currFaceID/2); 
            p = prepareFaceCalcs(p,localIdx,G,mcells(i,:),fracMatAreas(i),currFaceID, projectionCells);        
        elseif (distRight+tol)<distLeft
            %right
            currFaceID = 2;
            faceSelected(currFaceID) = 1;                        
            localIdx = 3*(i-1)+ceil(currFaceID/2);
            p = prepareFaceCalcs(p,localIdx,G,mcells(i,:),fracMatAreas(i),currFaceID, projectionCells);
        else           
            %2 Equidistant faces in the X-direction
            isEquiDist(1) = 1;
        end

        if (distFront+tol)<distBack
            %front
            currFaceID = 3;
            faceSelected(currFaceID) = 1;                        
            localIdx = 3*(i-1)+ceil(currFaceID/2);
            p = prepareFaceCalcs(p,localIdx,G,mcells(i,:),fracMatAreas(i),currFaceID, projectionCells);
        elseif (distBack+tol)<distFront
            %back
            currFaceID = 4;
            faceSelected(currFaceID) = 1;                        
            localIdx = 3*(i-1)+ceil(currFaceID/2);
            p = prepareFaceCalcs(p,localIdx,G,mcells(i,:),fracMatAreas(i),currFaceID, projectionCells);
        else           
            isEquiDist(2) = 1;
        end

        if (distBottom+tol)<distTop 
            %bottom
            currFaceID = 5;
            faceSelected(currFaceID) = 1;                        
            localIdx = 3*(i-1)+ceil(currFaceID/2);
            p = prepareFaceCalcs(p,localIdx,G,mcells(i,:),fracMatAreas(i),currFaceID, projectionCells);  
        elseif (distTop+tol)<distBottom
            %top
            currFaceID = 6;
            faceSelected(currFaceID) = 1;                        
            localIdx = 3*(i-1)+ceil(currFaceID/2);
            p = prepareFaceCalcs(p,localIdx,G,mcells(i,:),fracMatAreas(i),currFaceID, projectionCells);
        else           
            isEquiDist(3) = 1;
        end
        
        %Algorithm to choose continuous projection cells if 2 faces along
        %the same direction are equidistant from the fracture centroid
        lia = ismember(G.nnc.cells, mcells(i,:),'rows');
        nncIdx4MF = find(lia);
        planepoint = G.nnc.planepoint(nncIdx4MF,:);
        fracFaceNormal = G.nnc.normal(nncIdx4MF,:);
            
        numEquiDist = sum(isEquiDist);
        if numEquiDist==0
            %case 1: Frac centroid is closer to 3 faces than other 3 faces
            %No additional work to be done here.
            
        elseif numEquiDist==1
            %case 2: Frac centroid is closer to 2 faces than the opposite
            %2 faces. It is equidistant from two faces in 1 direction
            equiDistDir = find(isEquiDist); %face direction that is equidistant from frac centroid
%             selectedDir = find(~isEquiDist); %non-equidistant face directions from frac centroid
%             projCellIdx = find(projectionCell); %indices of the 2 selected projection cells
            projFaceIdx = find(faceSelected); %indices of the 2 selected projection faces
            projFaces = faceIdxs_i(projFaceIdx); %The 2 selected projection faces (global face numbers)
            
            nodeStartIdx = G.faces.nodePos(projFaces');   %indices of the first node of the 2 selected faces
            nodeEndIdx = G.faces.nodePos(projFaces'+1)-1; %indices of the last node of the 2 selected faces
            nodeIdxs_1 = G.faces.nodes(nodeStartIdx(1):nodeEndIdx(1)); %indices of nodes that bound 1st selected face
            nodeIdxs_2 = G.faces.nodes(nodeStartIdx(2):nodeEndIdx(2)); %indices of nodes that bound 2nd selected face
            commonNodeIdx_s = intersect(nodeIdxs_1,nodeIdxs_2); %return the 2 node indices common to the 2 selected faces
            commonNodeCoords = G.nodes.coords(commonNodeIdx_s,:);%return the XYZ coords of these two common nodes

            %find distance between the fracture plane and the 2 nodes common
            %to the 2 faces that were selected based on unequal distances along 2 axes
            normdist = pointplanedistance(commonNodeCoords,repmat(fracFaceNormal,2,1),repmat(planepoint,2,1));
            [~, idxMaxDist] = max(normdist);
%             nodeXYZ = commonNodeCoords(idxMaxDist,:);
               
            localIdxEqui = 3*(i-1)+equiDistDir;
            idx1to6 = 2*(equiDistDir-1)+idxMaxDist;
%             localIdxNonEqui = 3*(i-1)+selectedDir;
            if (idxMaxDist==1)
                %this means cell to the left, front or bottom is selected      
                projectionCell(equiDistDir) = projectionCells(idx1to6,1);         
                p.MMneighs(localIdxEqui,:) = [projectionCell(equiDistDir),mcells(i,1)];
            elseif (idxMaxDist==2)
                %this means cell to the right, back or top is selected
                projectionCell(equiDistDir) = projectionCells(idx1to6,2);         
                p.MMneighs(localIdxEqui,:) = [mcells(i,1), projectionCell(equiDistDir)];
            else
                printf('Error!, code should never get here. Wrong idxMaxDist for case 2 (numEquiDist==1) ');
            end
            
            p.MFcells(localIdxEqui,:) = [projectionCell(equiDistDir), mcells(i,2)];
            [p.Type{localIdxEqui}, p.Area(localIdxEqui), p.Trans(localIdxEqui), p.Dir(localIdxEqui)] = ...
                    computeNNCprojArea_n_Trans(G,p.MMneighs(localIdxEqui,:),p.MFcells(localIdxEqui,:),mcells(i,:),fracMatAreas(i));

            
        elseif numEquiDist==2
            %case 3: Frac centroid is closer to 1 face than the other opposite
            %face. It is equidistant from the other 4 faces in 2 directions
%             equiDistDir = find(isEquiDist); %face directions that are equidistant from frac centroid
            selectedDir = find(~isEquiDist); %non-equidistant face direction from frac centroid
%             projCellIdx = find(projectionCell); %index of the selected projection cell
%             selectedProjectionCell = projectionCell(projCellIdx);
            projFaceIdx = find(faceSelected); %index of the selected projection face
            projFace = faceIdxs_i(projFaceIdx); %The selected projection face (global face number)
            
            nodeStartIdx = G.faces.nodePos(projFace);   %index of the first node of the selected face
            nodeEndIdx = G.faces.nodePos(projFace+1)-1; %index of the last node of the selected face
            nodeIdxs = G.faces.nodes(nodeStartIdx:nodeEndIdx); %indices of nodes that bound selected face
            nodeCoords = G.nodes.coords(nodeIdxs,:);%return the XYZ coords of the selected face
            %find distance between all four nodes of the selected face and
            %fracture plane
            normdist = pointplanedistance(nodeCoords,repmat(fracFaceNormal,4,1),repmat(planepoint,4,1));
            
            [~, idxMaxDist] = max(normdist);
            
            %Note that nnc data for the selected face has already been
            %computed. We only need to focus on the other two sides that
            %meet the selected face.
            
            %compute nnc data for the other two sides that meet with the
            %selected left face at the selected node at index idxMaxDist
            if (selectedDir==1)              
                if (idxMaxDist==1)
                    %select front and bottom 
                    %front
                    currFaceID = 3;                        
                    localIdx = 3*(i-1)+ceil(currFaceID/2);
                    p = prepareFaceCalcs(p,localIdx,G,mcells(i,:),fracMatAreas(i),currFaceID, projectionCells);
                    %bottom
                    currFaceID = 5;                        
                    localIdx = 3*(i-1)+ceil(currFaceID/2);
                    p = prepareFaceCalcs(p,localIdx,G,mcells(i,:),fracMatAreas(i),currFaceID, projectionCells); 
                elseif (idxMaxDist==2)
                    %select back and bottom
                    %back
                    currFaceID = 4;                        
                    localIdx = 3*(i-1)+ceil(currFaceID/2);
                    p = prepareFaceCalcs(p,localIdx,G,mcells(i,:),fracMatAreas(i),currFaceID, projectionCells);
                    %bottom
                    currFaceID = 5;                        
                    localIdx = 3*(i-1)+ceil(currFaceID/2);
                    p = prepareFaceCalcs(p,localIdx,G,mcells(i,:),fracMatAreas(i),currFaceID, projectionCells);
                elseif (idxMaxDist==3)
                    %select back and top
                    %back
                    currFaceID = 4;                        
                    localIdx = 3*(i-1)+ceil(currFaceID/2);
                    p = prepareFaceCalcs(p,localIdx,G,mcells(i,:),fracMatAreas(i),currFaceID, projectionCells);
                    %top
                    currFaceID = 6;                        
                    localIdx = 3*(i-1)+ceil(currFaceID/2);
                    p = prepareFaceCalcs(p,localIdx,G,mcells(i,:),fracMatAreas(i),currFaceID, projectionCells);
                    
                elseif (idxMaxDist==4)
                    %select front and top
                    %front
                    currFaceID = 3;                        
                    localIdx = 3*(i-1)+ceil(currFaceID/2);
                    p = prepareFaceCalcs(p,localIdx,G,mcells(i,:),fracMatAreas(i),currFaceID, projectionCells);
                    %top
                    currFaceID =6;                        
                    localIdx = 3*(i-1)+ceil(currFaceID/2);
                    p = prepareFaceCalcs(p,localIdx,G,mcells(i,:),fracMatAreas(i),currFaceID, projectionCells); 
                else
                    printf('Error!, code should never get here. Wrong idxMaxDist for case 3 (numEquiDist==2) ');
                end
                
            elseif (selectedDir==2)
                if (idxMaxDist==1)
                    %select left and bottom 
                    %left
                    currFaceID = 1;                        
                    localIdx = 3*(i-1)+ceil(currFaceID/2);
                    p = prepareFaceCalcs(p,localIdx,G,mcells(i,:),fracMatAreas(i),currFaceID, projectionCells);
                    %bottom
                    currFaceID = 5;                        
                    localIdx = 3*(i-1)+ceil(currFaceID/2);
                    p = prepareFaceCalcs(p,localIdx,G,mcells(i,:),fracMatAreas(i),currFaceID, projectionCells);
                elseif (idxMaxDist==2)
                    %select left and top
                    %left 
                    currFaceID = 1;                        
                    localIdx = 3*(i-1)+ceil(currFaceID/2);
                    p = prepareFaceCalcs(p,localIdx,G,mcells(i,:),fracMatAreas(i),currFaceID, projectionCells);
                    %top
                    currFaceID = 6;                        
                    localIdx = 3*(i-1)+ceil(currFaceID/2);
                    p = prepareFaceCalcs(p,localIdx,G,mcells(i,:),fracMatAreas(i),currFaceID, projectionCells);
                elseif (idxMaxDist==3)
                    %select right and top
                    %right
                    currFaceID = 2;                        
                    localIdx = 3*(i-1)+ceil(currFaceID/2);
                    p = prepareFaceCalcs(p,localIdx,G,mcells(i,:),fracMatAreas(i),currFaceID, projectionCells);
                    %top
                    currFaceID = 6;                        
                    localIdx = 3*(i-1)+ceil(currFaceID/2);
                    p = prepareFaceCalcs(p,localIdx,G,mcells(i,:),fracMatAreas(i),currFaceID, projectionCells);
                elseif (idxMaxDist==4)
                    %select right and bottom
                    %right
                    currFaceID = 2;                        
                    localIdx = 3*(i-1)+ceil(currFaceID/2);
                    p = prepareFaceCalcs(p,localIdx,G,mcells(i,:),fracMatAreas(i),currFaceID, projectionCells);
                    %bottom
                    currFaceID = 5;                        
                    localIdx = 3*(i-1)+ceil(currFaceID/2);
                    p = prepareFaceCalcs(p,localIdx,G,mcells(i,:),fracMatAreas(i),currFaceID, projectionCells);
                else
                    printf('Error!, code should never get here. Wrong idxMaxDist for case 3 (numEquiDist==2) ');
                end
                
            elseif (selectedDir==3)
                if (idxMaxDist==1)
                    %select left and front
                    %left
                    currFaceID = 1;                        
                    localIdx = 3*(i-1)+ceil(currFaceID/2);
                    p = prepareFaceCalcs(p,localIdx,G,mcells(i,:),fracMatAreas(i),currFaceID, projectionCells);
                    %front
                    currFaceID = 3;                        
                    localIdx = 3*(i-1)+ceil(currFaceID/2);
                    p = prepareFaceCalcs(p,localIdx,G,mcells(i,:),fracMatAreas(i),currFaceID, projectionCells);                  
                elseif (idxMaxDist==2)
                    %select right and front
                    %right
                    currFaceID = 2;                        
                    localIdx = 3*(i-1)+ceil(currFaceID/2);
                    p = prepareFaceCalcs(p,localIdx,G,mcells(i,:),fracMatAreas(i),currFaceID, projectionCells);  
                    %front
                    currFaceID = 3;                        
                    localIdx = 3*(i-1)+ceil(currFaceID/2);
                    p = prepareFaceCalcs(p,localIdx,G,mcells(i,:),fracMatAreas(i),currFaceID, projectionCells);                   
                elseif (idxMaxDist==3)
                    %select right and back
                    %right
                    currFaceID = 2;                        
                    localIdx = 3*(i-1)+ceil(currFaceID/2);
                    p = prepareFaceCalcs(p,localIdx,G,mcells(i,:),fracMatAreas(i),currFaceID, projectionCells);
                    %back
                    currFaceID = 4;                        
                    localIdx = 3*(i-1)+ceil(currFaceID/2);
                    p = prepareFaceCalcs(p,localIdx,G,mcells(i,:),fracMatAreas(i),currFaceID, projectionCells);  
                elseif (idxMaxDist==4)
                    %select left and back
                    %left
                    currFaceID = 1;                        
                    localIdx = 3*(i-1)+ceil(currFaceID/2);
                    p = prepareFaceCalcs(p,localIdx,G,mcells(i,:),fracMatAreas(i),currFaceID, projectionCells);
                    %back
                    currFaceID = 4;                        
                    localIdx = 3*(i-1)+ceil(currFaceID/2);
                    p = prepareFaceCalcs(p,localIdx,G,mcells(i,:),fracMatAreas(i),currFaceID, projectionCells); 
                else
                    printf('Error!, code should never get here. Wrong idxMaxDist for case 3 (numEquiDist==2) ');
                end
                
            else
                printf('Code should never get here. Wrong selectedDir for case 3 (numEquiDist==2) ');
            end
            
        elseif numEquiDist==3
            %case 4: Frac centroid is equidistant from all 6 cell faces
            switch(lower(opt.defaultFace))
                case 'left'
                    isEquiDist(1) = 0; %randomly select left side as the selected side
                    faceSelected(1) = 1;
                    
                    selectedDir = find(~isEquiDist); %non-equidistant face direction from frac centroid
                    projFaceIdx = find(faceSelected); %index of the selected projection face
                    projFace = faceIdxs_i(projFaceIdx); %The selected projection face (global face number)

                    nodeStartIdx = G.faces.nodePos(projFace);   %index of the first node of the selected face
                    nodeEndIdx = G.faces.nodePos(projFace+1)-1; %index of the last node of the selected face
                    nodeIdxs = G.faces.nodes(nodeStartIdx:nodeEndIdx); %indices of nodes that bound selected face
                    nodeCoords = G.nodes.coords(nodeIdxs,:);%return the XYZ coords of the selected face
                    %find distance between all four nodes of the selected face and
                    %fracture plane
                    normdist = pointplanedistance(nodeCoords,repmat(fracFaceNormal,4,1),repmat(planepoint,4,1));
                    [~, idxMaxDist] = max(normdist);

                    %compute nnc data for the randomly selected left face
                    localIdxEqui = 3*(i-1)+selectedDir;     
                    projectionCell(selectedDir) = projectionCells(1,1);         
                    p.MMneighs(localIdxEqui,:) = [projectionCell(selectedDir),mcells(i,1)];

                    p.MFcells(localIdxEqui,:) = [projectionCell(selectedDir), mcells(i,2)];
                    [p.Type{localIdxEqui}, p.Area(localIdxEqui), p.Trans(localIdxEqui), p.Dir(localIdxEqui)] = ...
                            computeNNCprojArea_n_Trans(G,p.MMneighs(localIdxEqui,:),p.MFcells(localIdxEqui,:),mcells(i,:),fracMatAreas(i));

                    %compute nnc data for the other two sides that meet with the
                    %selected left face at the selected node at index idxMaxDist
                    if (idxMaxDist==1)
                        %select front and bottom 
                        %front
                        currFaceID = 3;                        
                        localIdx = 3*(i-1)+ceil(currFaceID/2);
                        p = prepareFaceCalcs(p,localIdx,G,mcells(i,:),fracMatAreas(i),currFaceID, projectionCells);
                        %bottom
                        currFaceID = 5;                        
                        localIdx = 3*(i-1)+ceil(currFaceID/2);
                        p = prepareFaceCalcs(p,localIdx,G,mcells(i,:),fracMatAreas(i),currFaceID, projectionCells); 
                    elseif (idxMaxDist==2)
                        %select back and bottom
                        %back
                        currFaceID = 4;                        
                        localIdx = 3*(i-1)+ceil(currFaceID/2);
                        p = prepareFaceCalcs(p,localIdx,G,mcells(i,:),fracMatAreas(i),currFaceID, projectionCells);
                        %bottom
                        currFaceID = 5;                        
                        localIdx = 3*(i-1)+ceil(currFaceID/2);
                        p = prepareFaceCalcs(p,localIdx,G,mcells(i,:),fracMatAreas(i),currFaceID, projectionCells);
                        
                    elseif (idxMaxDist==3)
                        %select back and top
                        %back
                        currFaceID = 4;                        
                        localIdx = 3*(i-1)+ceil(currFaceID/2);
                        p = prepareFaceCalcs(p,localIdx,G,mcells(i,:),fracMatAreas(i),currFaceID, projectionCells);
                        %top
                        currFaceID = 6;                        
                        localIdx = 3*(i-1)+ceil(currFaceID/2);
                        p = prepareFaceCalcs(p,localIdx,G,mcells(i,:),fracMatAreas(i),currFaceID, projectionCells);
                        
                    elseif (idxMaxDist==4)
                        %select front and top
                        %front
                        currFaceID = 3;                        
                        localIdx = 3*(i-1)+ceil(currFaceID/2);
                        p = prepareFaceCalcs(p,localIdx,G,mcells(i,:),fracMatAreas(i),currFaceID, projectionCells);
                        %top
                        currFaceID = 6;                        
                        localIdx = 3*(i-1)+ceil(currFaceID/2);
                        p = prepareFaceCalcs(p,localIdx,G,mcells(i,:),fracMatAreas(i),currFaceID, projectionCells);
                        
                    else
                        printf('Error!, code should never get here. Wrong idxMaxDist for case 4 (numEquiDist==3) ');
                    end
                case 'front'
%                     isEquiDist(2) = 0; %randomly select front side as the selected side
%                     faceSelected(3) = 1;

                    isEquiDist(2) = 0; %randomly select front side as the selected side
                    faceSelected(3) = 1;
                case 'bottom'
                    isEquiDist(3) = 0; %randomly select bottom side as the selected side
                    faceSelected(5) = 1;
                    
                    selectedDir = find(~isEquiDist); %non-equidistant face direction from frac centroid
                    projFaceIdx = find(faceSelected); %index of the selected projection face
                    projFace = faceIdxs_i(projFaceIdx); %The selected projection face (global face number)

                    nodeStartIdx = G.faces.nodePos(projFace);   %index of the first node of the selected face
                    nodeEndIdx = G.faces.nodePos(projFace+1)-1; %index of the last node of the selected face
                    nodeIdxs = G.faces.nodes(nodeStartIdx:nodeEndIdx); %indices of nodes that bound selected face
                    nodeCoords = G.nodes.coords(nodeIdxs,:);%return the XYZ coords of the selected face
                    %find distance between all four nodes of the selected face and
                    %fracture plane
                    normdist = pointplanedistance(nodeCoords,repmat(fracFaceNormal,4,1),repmat(planepoint,4,1));
                    [~, idxMaxDist] = max(normdist);

                    %compute nnc data for the randomly selected bottom face
                    localIdxEqui = 3*(i-1)+selectedDir;     
                    projectionCell(selectedDir) = projectionCells(5,1);         
                    p.MMneighs(localIdxEqui,:) = [projectionCell(selectedDir),mcells(i,1)];
                    p.MFcells(localIdxEqui,:) = [projectionCell(selectedDir), mcells(i,2)];
                    [p.Type{localIdxEqui}, p.Area(localIdxEqui), p.Trans(localIdxEqui), p.Dir(localIdxEqui)] = ...
                            computeNNCprojArea_n_Trans(G,p.MMneighs(localIdxEqui,:),p.MFcells(localIdxEqui,:),mcells(i,:),fracMatAreas(i));
         
                    %compute nnc data for the other two sides that meet with the
                    %selected left face at the selected node at index idxMaxDist
                    if (idxMaxDist==1)
                        %select left and front 
                        %left  ( %left/right => +1; front/back => +2; bottom/top => +3;)
                        currFaceID = 1;                        
                        localIdx = 3*(i-1)+ceil(currFaceID/2);
                        p = prepareFaceCalcs(p,localIdx,G,mcells(i,:),fracMatAreas(i),currFaceID, projectionCells);
                        %front
                        currFaceID = 3;                        
                        localIdx = 3*(i-1)+ceil(currFaceID/2);
                        p = prepareFaceCalcs(p,localIdx,G,mcells(i,:),fracMatAreas(i),currFaceID, projectionCells);
                        
                    elseif (idxMaxDist==2)
                        %select right and front
                        %right
                        currFaceID = 2;                        
                        localIdx = 3*(i-1)+ceil(currFaceID/2);
                        p = prepareFaceCalcs(p,localIdx,G,mcells(i,:),fracMatAreas(i),currFaceID, projectionCells);
                        %front
                        currFaceID = 3;                        
                        localIdx = 3*(i-1)+ceil(currFaceID/2);
                        p = prepareFaceCalcs(p,localIdx,G,mcells(i,:),fracMatAreas(i),currFaceID, projectionCells);
                        
                    elseif (idxMaxDist==3)
                        %select right and back
                        %right
                        currFaceID = 2;                        
                        localIdx = 3*(i-1)+ceil(currFaceID/2);
                        p = prepareFaceCalcs(p,localIdx,G,mcells(i,:),fracMatAreas(i),currFaceID, projectionCells);
                        %back
                        currFaceID = 4;                        
                        localIdx = 3*(i-1)+ceil(currFaceID/2);
                        p = prepareFaceCalcs(p,localIdx,G,mcells(i,:),fracMatAreas(i),currFaceID, projectionCells);
                        
                    elseif (idxMaxDist==4)
                        %select left and back
                        %left
                        currFaceID = 1;                        
                        localIdx = 3*(i-1)+ceil(currFaceID/2);
                        p = prepareFaceCalcs(p,localIdx,G,mcells(i,:),fracMatAreas(i),currFaceID, projectionCells);
                        %back
                        currFaceID = 4;                        
                        localIdx = 3*(i-1)+ceil(currFaceID/2);
                        p = prepareFaceCalcs(p,localIdx,G,mcells(i,:),fracMatAreas(i),currFaceID, projectionCells);

                    else
                        printf('Error!, code should never get here. Wrong idxMaxDist for case 4 (numEquiDist==3) ');
                    end
                otherwise
                    error('Unknown case')
            end

        else
            %code should never get here
            printf('Code should never get here! numEquiDist is wrong!')
        end


    end
    
    toc
    
    %Modification to remove projection matrices on domain boundaries
    idx2keep = all(p.MMneighs~=0, 2);
    p.MMneighs = p.MMneighs(idx2keep,:);
    p.MFcells = p.MFcells(idx2keep,:);
    p.Type = p.Type(idx2keep,:);
    p.Area = p.Area(idx2keep,:);
    p.Dir = p.Dir(idx2keep,:);
    
    p.Trans = p.Trans(idx2keep,:);
    G.nnc.T = [G.nnc.T; p.Trans];
    
    %% Compute union area
    %Code to extract the nodes of a list of cells (lines 68-74 in EDFMgrid.m)
    if ~isempty(mcells)
        [cn,cpos]=gridCellNodes(G,mcells(1,2):mcells(end,2));
    end
    numFcells = numel(mcells(:,2)); % mcells {matrix_cell,frac_cell}
    mcellnodes_par=cell(numFcells,1);
    for i=1:numFcells
        mcellnodeind=cn(cpos(i):(cpos(i+1)-1)); %idx of all nodes of each cell
        mcellnodes=G.nodes.coords(mcellnodeind,:);
        mcellnodes_par{i}=mcellnodes;
    end
    
% https://www.mathworks.com/matlabcentral/answers/336500-finding-the-indices-of-duplicate-values-in-one-array
    [pMMneighsUniq, uniqIdx, ic2] = unique(p.MMneighs,'rows','stable');
    duplicateMatIdx2 = setdiff( 1:numel(p.MMneighs(:,1)), uniqIdx); %This gives you the indices of the 2nd and subsequent occurrences of a matrix cell
    p.Area = p.Area(uniqIdx); %these are the areas for the unique pairs of pM-M cells 
    repIdx = ic2(duplicateMatIdx2);
    
    pMMrepeated = p.MMneighs(duplicateMatIdx2,:);
    pDirrepeated = p.Dir(duplicateMatIdx2);
    [~, ia3, ~] = unique(pMMrepeated,'rows','stable');
    
    %Loop over each matrix cell with multiple fractures
    for mi=1:numel(ia3)
        pMMID = pMMrepeated(ia3(mi),:);
        pDir = pDirrepeated(ia3(mi));
        lia2 = ismember(p.MMneighs,pMMID,'rows');
        fracIDs = p.MFcells(lia2,2);
        lia = ismember(mcells(:,2),fracIDs);
        FracPolyVerts=mcellnodes_par(lia);
        repIdx2 = repIdx(ia3(mi));
        
        area =getUnionProjectArea(FracPolyVerts,pDir,tol); %Proj Area in X,Y,Z
        if isnan(area)
            area =getUnionProjectArea(FracPolyVerts,pDir);
        else
            p.Area(repIdx2) = area;
        end
        
    end


    p.MMneighs = pMMneighsUniq; %p.MMneighs(uniqIdx,:); %these are same
    p.MFcells = p.MFcells(uniqIdx,:);
    p.Type = p.Type(uniqIdx,:);
    p.Dir = p.Dir(uniqIdx,:);
    G.nnc.pDir = p.Dir;
    
    G.nnc.cells = [G.nnc.cells; p.MFcells];
    G.nnc.type = [G.nnc.type; p.Type];
    G.nnc.area = [G.nnc.area; p.Area];
    G.nnc.pMMneighs = p.MMneighs; %This is technically not an NNC but is 
    %saved as such for convenience. It will be used to compute transmultpEDFM
    % for modifying existing projected Matrix - Matrix Transmissibilities
    
 
end

function [pType, pArea, pTrans, projDir] = computeNNCprojArea_n_Trans(G,pMMneighs,pMFcells,mcells,fracMatArea)
    if (pMFcells(1)~=0)
        %if matrix host cell is on the boundary, projection is
        %unncessary. Initialize all nnc arrays correctly in this case
        pType = 'fracpmat';
        lia = ismember(G.faces.neighbors,pMMneighs,'rows');
        globalIdx4pMM = find(lia); 
        lia = ismember(G.nnc.cells, mcells,'rows');
        nncIdx4MF = find(lia);  
        
        pArea = abs(fracMatArea*(dot(G.faces.normals(globalIdx4pMM,:),G.nnc.normal(nncIdx4MF,:)))...
            /(norm(G.faces.normals(globalIdx4pMM,:))*norm(G.nnc.normal(nncIdx4MF,:))));
        % Projection cell is pMFcells(1) 
        % current fraccellcentroid = G.cells.centroids(mcells(2),:)
        dist_pMF = norm(G.cells.centroids(pMFcells(1),:)-G.cells.centroids(mcells(2),:));
%       K_MpF = 2.0*G.rock.perm(pMFcells(localIdx,1))*G.rock.perm(pMFcells(localIdx,2))/(G.rock.perm(pMFcells(localIdx,1))+G.rock.perm(pMFcells(localIdx,2)));
%        K_pMF = harmmean(G.rock.perm(pMFcells));
        K_pMF = length(pMFcells)/sum(1./G.rock.perm(pMFcells)); %is this right?
        pTrans = K_pMF*pArea/dist_pMF; %pTrans is T_pM-F
        projDir = find(G.faces.normals(globalIdx4pMM,:));
    else
        %it means pM cell is on the boundary of the domain, and we don't
        %need any projection
        pType = 'ProjMatrix outside';
        pArea = 0;
        pTrans = 0;    
        projDir = 0;
    end 

end

function [p] = prepareFaceCalcs(p,localIdx,G,mcells,fracMatArea,currFaceID,projectionCells)                  
    if (mod(currFaceID,2)==1)
        %this means cell to the left, front or bottom is selected
        projectionCell = projectionCells(currFaceID, 1);
        p.MFcells(localIdx,:) = [projectionCell, mcells(2)];
        p.MMneighs(localIdx,:) = [projectionCell, mcells(1)];
    else
        %this means cell to the right, back or top is selected
        projectionCell = projectionCells(currFaceID, 2);
        p.MFcells(localIdx,:) = [projectionCell, mcells(2)];
        p.MMneighs(localIdx,:) = [mcells(1), projectionCell];
    end
    [p.Type{localIdx}, p.Area(localIdx), p.Trans(localIdx), p.Dir(localIdx)] = ...
            computeNNCprojArea_n_Trans(G,p.MMneighs(localIdx,:),p.MFcells(localIdx,:),mcells,fracMatArea);                 
end

function [areas]=getUnionProjectArea(FracPolyVerts,pDir,tol)
    % Compute union projection area with the multiple fractures in x,y,z dirs

    NumFracs=numel(FracPolyVerts);

    %Convert box frac to plane frac OMO:will this work for any orentatn?
    FracNorms=zeros(NumFracs,3);
    for fi=1:NumFracs
        poly_box=FracPolyVerts{fi};
        numVert=size(poly_box,1);
        poly_plane=(poly_box(1:numVert/2,:)+poly_box(numVert/2+1:end,:))/2.0;
        FracPolyVerts{fi}=poly3Dccw(poly_plane);
        FracNorms(fi,:)=poly3Dnormal(FracPolyVerts{fi});
    end

    %Frac polygon union
    polyvec_x=[];
    polyvec_y=[];
    polyvec_z=[];

    for fi=1:NumFracs
        X=FracPolyVerts{fi}(:,1);
        Y=FracPolyVerts{fi}(:,2);
        Z=FracPolyVerts{fi}(:,3);

        if abs(FracNorms(fi,1))>tol %Not prependicular to YZ Plane
           polyvec_x=[polyvec_x polyshape(Y,Z)]; %YZ Plane Projection
        end
        if abs(FracNorms(fi,2))>tol %Equal to FracNorms(fi,2)~=0 
           polyvec_y=[polyvec_y polyshape(X,Z)]; %XZ Plane Projection
        end
        if abs(FracNorms(fi,3))>tol  
           polyvec_z=[polyvec_z polyshape(X,Y)]; %XY Plane Projection
        end
    end

    %Compute the union areas
%     areas=[0.0,0.0,0.0];
    areas=0.0;
    if ~isempty(polyvec_x) && pDir==1
        polyout_x = union(polyvec_x);
        areas = area(polyout_x); %Equal to sum(area(regions(polyout_x)))
    end

    if ~isempty(polyvec_y) && pDir==2
        polyout_y = union(polyvec_y);
        areas = area(polyout_y);
    end

    if ~isempty(polyvec_z) && pDir==3
        polyout_z = union(polyvec_z);
        areas = area(polyout_z);
    end


    %debug
    %figure
    %plot(polyvec_x)

    %figure
    %plot(polyout_x)

    %close all
    %plot3(xyz(:,1),xyz(:,2),xyz(:,3))

end

function [norm_vec] = poly3Dnormal(xyz)
%Compute unit normal vector for each fracture plane
% Plane normals
points = [xyz;xyz(1,:)];
diffp = diff(points,1);
normal = cross(diffp(1,:), diffp(2,:));
norm_vec = normal/norm(normal);
end

function [xyz_order] = poly3Dccw(xyz)
%https://www.mathworks.com/matlabcentral/answers/429265-how-to-
%order-vertices-of-a-flat-convex-polygon-in-3d-space-along-the-edge
    xyzc = mean(xyz,1);
    P = xyz - xyzc;
    [~,~,V] = svd(P,0);
    [~,is] = sort(atan2(P*V(:,1),P*V(:,2)));
    xyz_order = xyz(is([1:end]),:);
end
