function altRotate(varargin)
%Alternative function replacing rotate3d for rotation of 3d figures
%
% DESCRIPTION:
%  This function is called for rotation in some mrst gui plots. The
%  function is intended as a drop-in replacement of rotate3d, which is
%  normally called during object rotation (patch plots and so on). The
%  reason for the replacement is that the default rotation is very slow,
%  due to calls to drawnow as well as several calls to check for legacy
%  functionality.
%
%  altRotate is automatically used by functions where the primary interest
%  is visualization of a 3D reservoir model where rotate3d is too slow to
%  be usable at all.
%
%  To replace the rotate button in an existing plot, see "fastRotateButton"

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


persistent oldMode;

if nargin == 0
    h = gcf;
else
    h = varargin{1};
end

isHG1 = isnumeric(h);

h = double(h);

if nargin > 1
    state = varargin{2};
else
    state = true;
end


if isnan(h) || isempty(h)
    h = gcf;
end

if state
    set(h, 'WindowButtonDownFcn', @WindowButtonDownFcn);
    set(h, 'WindowButtonUpFcn', @WindowButtonUpFcn);
else
    set(h, 'WindowButtonDownFcn', []);
    set(h, 'WindowButtonUpFcn', []);
end


function WindowButtonDownFcn(src, event)
    if isHG1
        oldMode = get(gca, 'DrawMode');
        set(gca, 'DrawMode', 'fast');
    else
        oldMode = get(gca, 'SortMethod');
        set(gca, 'SortMethod', 'childorder');
    end
    pt = get(h, 'CurrentPoint');
    pt = pt(1,:);

    [az, el] = view();

    % Cap to 90 degrees elevation in either direction to avoid strange axis
    % effects
    rotatefn = @(a)  min(max(...
                     ([az, el] - (a(1,:) - pt) + [0 90]),...
                     [-inf, 0]), [inf, 180]) - [0 90];
    ax = gca;


    set(h, 'WindowButtonMotionFcn',...
           @(varargin) set(ax, 'View', rotatefn(get(h, 'CurrentPoint'))));

end

function WindowButtonUpFcn(src, event)
    set(h, 'WindowButtonMotionFcn', []);
    if isHG1
        set(gca, 'DrawMode', oldMode);
    else
        set(gca, 'SortMethod', oldMode);
    end
    drawnow
end


end

