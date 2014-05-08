function [xyz] = findMassCenter(g, resSol)
%{
#COPYRIGHT#
%}

xyz = sum(bsxfun(@times, resSol.s(:,1), g.cells.centroids))/sum(resSol.s(:,1));


end
