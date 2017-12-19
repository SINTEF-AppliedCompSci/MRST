function varargout=plotEDFMgrid(G)

figure; 
plotGrid(G,1:G.Matrix.cells.num,'facealpha',0);
plotGrid(G,(G.Matrix.cells.num+1):G.cells.num,'facealpha',0.75,'Facecolor','y');  

view(15,20);
axis equal tight;





end