function x = discspace1(d1, d2, N, refine, type)
%Discretize 1D space with options for refinement
%
% SYNOPSIS:
%   x = discspace1(start, stop, no_pts)
%   x = discspace1(start, stop, no_pts, refine)
%   x = discspace1(start, stop, no_pts, refine, type)
%
% REQUIRED PARAMETERS:
%   d1 - Start of interval.
%
%   d2 - End of interval. Function supports d1 > d2 if refine output is
%        modified accordingly.
%
%   N  - Number of points to discretize distance between d1 and d2. Note
%        that if the "constant" refinement is used, this number used as a
%        baseline for refinement, so that the actual number of points will
%        be larger. For other variants, the number of outputs is exact.
%
% OPTIONAL PARAMETERS:
%
%   refine - A N by 2 matrix that contains the start of refinement
%            intervals in the first column and a refinement factor in the
%            second. For example, for the interval [0, 1], refine could be:
%            [0.2,  4; % Use factor 1 from 0 to 0.2 (implicit) and 4 until
%             0.5,  3; % 0.5, where 3 is used as the factor until the next
%             0.75, 1] % point 0.75, where the default density is used.
%
%  type    - Either constant, where each interval is partitioned
%            separately, or a valid interpolation type for interp1.
%            Default: Linear.
%
% RETURNS:
%   x      - Partitioned interval.
% EXAMPLE:
%     figure; hold on
%     refine = [(1/5)*2*pi, 1; ...
%               (2/5)*2*pi, 5; ...
%               (3/5)*2*pi, 1]; % Refine interval of 2/3*2*pi -> 3/5 *2*pi
%     x = discspace1(0, 2*pi, 20, refine);
%     plot(x, sin(x), '.')
%     x = discspace1(0, 2*pi, 20, refine, 'pchip');
%     plot(x, sin(x), 'o')
%
% NOTE:
%   With three inputs, this function is identical to linspace.
%
% SEE ALSO:
%   linspace

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

    if nargin < 5
        type = 'linear';
        if nargin < 4
            refine = [];
        end
    end
    if d2 < d1
        refine = refine(end:-1:1, :);
        x = discspace1(d2, d1, N, refine, type);
        x = x(end:-1:1);
        return
    end
    assert(d1 ~= d2);
    if isempty(refine)
        x = linspace(d1, d2, N);
    else
        assert(all(diff(refine(:, 1)) > 0), 'First column of refinement must be monotone.');
        assert(size(refine, 2) == 2, 'Refinement must be N by 2.');
        assert(refine(end, 1) < d2, 'Refinement must end before final point');
        if refine(1, 1) > d1
            refine = [d1, 1; refine];
        end
        n = size(refine, 1);
        switch lower(type)
            case 'constant'
                Dx = d2 - d1;
                xc = cell(1, n);
                for i = 1:n
                    % Figure out roughly how many points interval would
                    % need if linspaced
                    if i == n
                        next = d2;
                    else
                        next = refine(i+1, 1);
                    end
                    current = refine(i, 1);
                    factor = refine(i, 2);
                    dx = next - current;
                    num = ceil(factor*N*dx/Dx);
                    
                    delta = dx/num;
                    xc{i} = linspace(current, next - delta, num);
                end
                x = horzcat(xc{:}, d2);
            otherwise
                X = refine(:, 1);
                Xs = (X - d1)./(d2 - d1);
                Xs = [Xs; 1];
                delta = diff([X; d2]);
                Y = cumsum([0; refine(:, 2).*delta]);
                Y = Y./Y(end);
                
                new = linspace(Y(1), Y(end), N);
                try
                    xs = interp1(Y, Xs, new, type);
                catch
                    error('Unknown refinement ''%s''. Valid options: constant, linear, spline.', type);
                end
                x = (d2 - d1)*xs + d1;
        end
    end
end
