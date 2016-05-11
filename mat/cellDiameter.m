function hK = cellDiameter(X)
%--------------------------------------------------------------------------
%   Calculates the longest distance between any two points in X.
%
%   SYNOPSIS:
%       hK = cellDiameter(X)
%
%   REQUIRED PARAMETERS:
%       X   - n x d matrix of points, where X are points in R^d.
%
%   RETURNS:
%       hK  - Cell diameter
%-----------------------------------------------------------------ØSK-2016-

%{
   Copyright (C) 2016 Øystein Strengehagen Klemetsdal. See Copyright.txt
   for details.
%}

n = size(X,1);
hK = 0;
for i = 1:n
    for j = i+1:n
        hK = max(norm(X(i,:)-X(j,:),2),hK);
    end
end
    
end