function G= reduceEDFMgridsize(G)
%REDUCEEDFMGRIDSIZE Reduces size of EDFM grid
%   Saves only cell/face/node numbers and global start indices in G.Matrix 
%   and G.FracGrid. Note that at the moment, this is irreversible since no
%   function has been written to reverse this process.

% G.matrixcellsnum=G.Matrix.cells.num;
% G=rmfield(G,{'Matrix','FracGrid'});

newMatrix=struct;
newMatrix.cells.num=G.Matrix.cells.num;
newMatrix.faces.num=G.Matrix.faces.num;
newMatrix.nodes.num=G.Matrix.nodes.num;
newMatrix.physdim=G.Matrix.physdim;
newMatrix.cartDims=G.Matrix.cartDims;

for i=1:numel(fieldnames(G.FracGrid))
    fieldname=['Frac',num2str(i)];
    
    newFracGrid.(fieldname).points=G.FracGrid.(fieldname).points;
    newFracGrid.(fieldname).cells.num=G.FracGrid.(fieldname).cells.num;
    newFracGrid.(fieldname).cells.start=G.FracGrid.(fieldname).cells.start;
    newFracGrid.(fieldname).faces.num=G.FracGrid.(fieldname).faces.num;
    newFracGrid.(fieldname).faces.start=G.FracGrid.(fieldname).faces.start;
    newFracGrid.(fieldname).nodes.num=G.FracGrid.(fieldname).nodes.num;
    newFracGrid.(fieldname).nodes.start=G.FracGrid.(fieldname).nodes.start;
    newFracGrid.(fieldname).matrixnnc=G.FracGrid.(fieldname).matrixnnc;
    newFracGrid.(fieldname).fracgridnnc=G.FracGrid.(fieldname).fracgridnnc;
end

G.Matrix=newMatrix;
G.FracGrid=newFracGrid;


end

