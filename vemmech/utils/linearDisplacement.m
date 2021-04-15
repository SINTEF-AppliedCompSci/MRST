function [disp, names] = linearDisplacement(dim)
%
%
% SYNOPSIS:
%   function [disp, names] = linearDisplacement(dim)
%
% DESCRIPTION:
%
% PARAMETERS:
%   dim - Dimension 
%
% RETURNS:
%   disp  - 
%   names - 
%
% EXAMPLE:
%
% SEE ALSO:
%

%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

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

% Compute the constant solutions
    disp = {}; names = {};
    for i = 1 : dim
        const = zeros(1, dim);
        const(i) = 1;
        disp{end + 1} = @(cc) repmat(const, size(cc, 1), 1);
        names{end + 1} = ['cons_disp_', num2str(i)];
    end

    % Compute the linear solutions
    for i = 1 : dim
        const = zeros(1, dim);
        const(i) = 1;
        disp{end + 1} = @(cc) repmat(bsxfun(@times, cc, const), 1);
        names{end + 1} = ['lin_disp_', num2str(i)];
    end

    if (dim == 2)
        disp{end + 1} = @(cc) cc(:, [2, 1]);
        names{end + 1} = 'steatch_xy';
        disp{end + 1} = @(cc) repmat(bsxfun(@times, cc(:, [2, 1]), [1, -1]), 1);
        names{end + 1} = 'rot_xy';
    else
        assert(dim == 3)
        % make permutations
        perm = [2 1 3;...
                1 3 2;...
                3 2 1]
        lnames = {'xy', 'yz', 'xz'};
        rigid = [1 -1 0;...
                 0  1 -1;...
                 1   0 -1];
        for i = 1 : 3
            % rigid rotation
            disp{end + 1} = @(cc) repmat(bsxfun(@times, cc(:, perm(i, :)), rigid(i, :)), 1);
            names{end + 1} = ['rot_', lnames{i}];
            % streatching
            disp{end + 1} = @(cc) repmat(bsxfun(@times, cc(:, perm(i, :)), abs(rigid(i, :))), 1);
            names{end + 1} = ['streatch_', lnames{i}];
            % streatching
        end    
    end
end