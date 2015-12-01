function [ Xcoord, Ycoord ] = getXYcoords(Gt, cellIndex)

    Xcoord = Gt.cells.centroids(cellIndex,1);
    Ycoord = Gt.cells.centroids(cellIndex,2);

end
