function fastRotateButton(varargin)
%Hook in faster rotation for 3D plots.
%
% SYNOPSIS:
%   fastRotateButton();
%
% DESCRIPTION:
%   This function replaces the rotate3d function with a much faster local
%   version specialized for 3d visualization. Used internally by mrst and
%   is subject to change.
%
% REQUIRED PARAMETERS:
%   None.
%
% RETURNS:
%   Nothing.
%
% NOTE:
%   fastRotateButton is primarily intended for internal use. It can be used
%   quite generally to speed up rotation, but no guarantees are made with
%   regards to parity with the original rotate3d function.
%
%
% EXAMPLE:
%    G = cartGrid([100, 100, 1]);
%    plotGrid(G);
%    fastRotateButton();
%
% SEE ALSO:
%   `altRotate`, `plotToolbar`

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

    if nargin > 0
        h = varargin{1};
    else
        h = gcf();
    end

    ht = findall(h,'tag','FigureToolBar');
    hRotate = findall(ht, 'tag', 'Exploration.Rotate');
    set(hRotate, 'ClickedCallback', @(varargin) rotateFast(varargin{:}));
end

function rotateFast(src, event, varargin)
    f = get(get(src, 'parent'), 'parent');
    activateuimode(f, '');
    if strcmpi(get(src, 'State'), 'on')
        altRotate(f, true);
        setptr(f,'rotate');
    else
        altRotate(f, false);
        setptr(f,'arrow');
    end
end
