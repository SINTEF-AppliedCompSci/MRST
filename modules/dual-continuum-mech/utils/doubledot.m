function result = doubledot(A, B, G)
%
% SYNOPSIS:
%   result = doubledot(A,B)
%
% DESCRIPTION:
%   Function to perform double dot product (or tensor contraction) on tensors
%   A and B, where tensors have been written in the oneline Voigt notation 
%   form used in the mechanics scripting. 
%
%   Note, this function relies on the fact that the tensors are symmetric
%
% PARAMETERS:
%   A, B - Tensors of rank 2 or 4
%   G - Grid structure
%
% RETURNS:
%   result 
%
% EXAMPLE:
%
% SEE ALSO:
%
%{
Copyright 2009-2020 SINTEF ICT, Applied Mathematics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}
if G.griddim == 3
    nlin = 6;
else
    nlin = 3;
end

% Check which tensor is of rank 4 (if not both)
    if size(A,2) == nlin^2 && size(B,2) == nlin^2
        rank4_A = A;
        rank4_B = B;
    elseif size(A,2) == nlin^2 
        rank4 = A;
        rank2 = B;  
    elseif size(B,2) == nlin^2
        rank4 = B;
        rank2 = A;
    else
        % must both be of rank 2, contraction then produces a scalar
        result = sum(A.*B, 2);
        return
    end

    try 
        rank4_A = reshape(permute(reshape(rank4_A, [], nlin), [2,1]), [], nlin);
        rank4_B = kron(rank4_B, ones(nlin,1));
        result = rank4_B .* repmat(rank4_A, 1, nlin);
        result = sum(reshape(permute(reshape(result, [], nlin),...
                                    [2,1]),...
                                    [], nlin),...
                                    2);
        result = reshape(result, nlin^2,[])';
    catch
        result = rank4 .* repmat(rank2, 1, nlin);
        result = sum(reshape(permute(reshape(result, [], nlin),...
                                    [2,1]),...
                                    [], nlin),...
                                    2);
        result = reshape(result, nlin,[])';
    end
    
end