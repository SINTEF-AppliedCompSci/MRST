function Gplot = createGplot(G,aplot,fracplanes)

Gplot = G;

cellstarts = zeros(length(fields(G.FracGrid)),1);
for i=1:length(fields(G.FracGrid))
    fieldname = ['Frac',num2str(i)];
    cellstarts(i) = G.FracGrid.(fieldname).cells.start;
end

for i=G.Matrix.cells.num+1:G.cells.num
    fgridnum = find(i>=cellstarts,1,'last');
    fieldname = ['Frac',num2str(fgridnum)];
    fracplanenum = G.FracGrid.(fieldname).fracplanenumber;
    normal = fracplanes(fracplanenum).normal;
    
    centroid = G.cells.centroids(i,:);
    
    cn = gridCellNodes(G,i);
    ncoords = G.nodes.coords(cn,:);
    
%     for j=1:length(cn)
%         perpdist = dot(ncoords(j,:)-centroid,normal);
%         ncoords(j,:) = ncoords(j,:) + (aplot/2 - perpdist)*normal;
%     end
    
    perpdist=dot((ncoords-centroid)',repmat(normal,length(cn),1)')';
    ncoords = ncoords - perpdist.*normal + perpdist./abs(perpdist)*aplot/2.*normal;
    
    Gplot.nodes.coords(cn,:) = ncoords;    
end

end