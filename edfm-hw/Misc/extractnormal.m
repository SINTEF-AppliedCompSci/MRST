function fracnormal=extractnormal(G,fracplanes,cellnum)
% Helper function to extract normal direction of fracture cell with global
% cell number 'cellnum'

fgnum=findFracGridnum(G,cellnum);

fracnormal=fracplanes(G.FracGrid.(['Frac',num2str(fgnum)]).fracplanenumber).normal;

end