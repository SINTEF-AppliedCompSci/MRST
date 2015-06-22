function [fraCells, remove] = markcells_old(G,fracplanes)
remove = [];
fraCells = cell(numel(fracplanes),1);
Sign = establishSign(G, fracplanes);
[cnodes, cnmap] = gridCellNodes(G,1:G.cells.num);
%
[fnodes,fnmap] = gridFaceNodes(G,1:G.faces.num);
nodesPerFace = diff(fnmap);
%
[cfaces,cfmap] = gridCellFaces(G,1:G.cells.num);
facesPerCell = diff(cfmap);
%
if any(diff(nodesPerFace)) % Hybrid/Composite Grid
    faceNodes = zeros(G.faces.num,max(nodesPerFace));
    cellFaces = zeros(G.cells.num,max(facesPerCell));
    for i = 1:G.faces.num
        temp = fnodes(fnmap(i):fnmap(i+1)-1).';
        faceNodes(i,1:numel(temp)) = temp;
    end
    for i = 1:G.cells.num
        temp = cfaces(cfmap(i):cfmap(i+1)-1).';
        cellFaces(i,1:numel(temp)) = temp;
    end
else
    faceNodes = reshape(fnodes,nodesPerFace(1),[]).';
    cellFaces = reshape(cfaces,facesPerCell(1),[]).';
end
count = 0;
for i = 1:numel(fracplanes)
    count = count+1;
    points = fracplanes(i).points;
    inside = sum(Sign(i).NodeSign(:,2:end),2)==size(points,1); % Nodes inside bounding planes
    bsum = Sign(i).NodeSign(:,2:end)>=0; % Points lying on boundary planes
    onplane = Sign(i).NodeSign(:,1) == 0; % Points on plane
    nodesOfInterest  = find(sum([inside,onplane,bsum],2) >= size(points,1)+1);
    if isempty(nodesOfInterest)
        warning(['Fracture ',num2str(i),' is not a long fracture and ',...
            'will be removed from further calculations.']);
        remove = [remove;i]; %#ok
        count = count-1;
        continue
    end
    isFace = ismember(faceNodes,nodesOfInterest);
    facesOfInterest = find(sum(isFace,2)>0); % faces with atleast 1 nodeOfInterest
    isCell = ismember(cellFaces,facesOfInterest);
    cellsOfInterest = find(sum(isCell,2)>0); % cells with atleast 1 faceOfInterest
    temp = [];
    for j = 1:numel(cellsOfInterest)
        ind = cellsOfInterest(j);
        nind = cnodes(cnmap(ind):cnmap(ind+1)-1);
        if numel(unique(Sign(i).NodeSign(nind,1)))>=2
            foi = cfaces(cfmap(ind):cfmap(ind+1)-1);
            isit = checkFracPenetration(G,fracplanes(i),ind,Sign(i),...
                foi,faceNodes,fnodes,fnmap);
            if isit
                temp = [temp;ind]; %#ok
            end
        end
    end 
    fraCells{count,1} = unique(temp);
end

return


    
            