function T_out = adjustT(T_in, cellpairs)
% Adjust T_in for surface area. Based on cellpairs, if a fracture cell
% intersects only one matrix cell, then the transmissibility between the
% two cells will be doubled since 'surface area = 2 x intersection area'.
% If a fracture cell intersects two matrix cells, then the transmissibility
% remains unchanged. Take note that only one pass of this function should
% be used. Refer to Tene et al. (2017) pEDFM for information on surface
% area.
%
% SYNOPSIS:
%   T_out = adjustT(T_in, cellpairs)
%
% REQUIRED PARAMETERS:
%   T_in    - Nx1 array of transmissibilities corresponding to each row in
%             cellpairs. T_in(i) is the transmissibility between
%             cellpairs(i,:) prior to adjustment for surface area.
%
%   cellpairs - Nx2 array with each row containing the cell indices of two
%               intersecting cells. cellpairs(i,1) corresponds to a matrix
%               cell and cellpairs(i,2) corresponds to a fracture cell.
%
% RETURNS:
%   T_out - New Nx1 array of transmissibilities adjusted for surface area

% Find indices in cellpairs(:,2) that are unique (non-repeated). These
% indices correspond to fracture cells contained within matrix cells.
bin = accumarray(cellpairs(:,2),1);
assert(all(bin<=2),'Found fracture cells intersecting more than 2 matrix cells. Check for bug.');
ind = ismember(cellpairs(:,2),find(bin==1));

% Multiply transmissibilities corresponding to those indices by 2
T_in(ind) = T_in(ind)*2;

T_out = T_in;


end