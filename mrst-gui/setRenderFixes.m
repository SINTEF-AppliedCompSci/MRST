function setRenderFixes(varargin)
%Set fixes which may reduce crashes in MATLAB plots.
%
% SYNOPSIS:
%   setRenderFixes()
%   setRenderFixes(gcf)
%   setRenderFixes(true)
%   setRenderFixes(gcf, false)
%
% DESCRIPTION:
%   This utility sets two figure and axes properties to values which are
%   known to remove some crashes and oddball behavior with MATLAB plots. If
%   you are not experiencing problems, this file should *NOT* be used.
%
%   The problems it attempts to fix are:
%          - Issued with multi monitor setups when figures become
%          unresponsive accompanied by a stack trace with references to jwt
%          errors. The fix is to set default renderer to opengl, which
%          enables plots to appear on any screen without strange behavior.
%          It should be noted that OpenGL plots have several limitations,
%          for instance missing support for logarithmic axes.
%
%          - Hard crash in MATLAB when plotting more complex figures
%          involving transparent patches, such as large reservoir cases.
%          The solution is to set drawmode to fast on the axes, but this
%          can have strange side effects when rotating objects later on.
%
%   Fixes can be applied and removed, either globally or per figure basis.
%   See next section.
%
% OPTIONAL PARAMETERS:
%
%   h      - Handle to figure fixes should be applied to. Defaults to 0
%            (root figure, defaults are applied system wide (!!!) ).
%
%   status - Boolean. If true, fixes are applied (Default). If false, fixes
%            are removed.
%
% RETURNS:
%     Nothing. Changes state of one or more figures.

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


    h = 0;
    status = true;

    for i = 1:numel(varargin)
        v = varargin{i};
        if islogical(v)
            status = v;
        elseif ishandle(v)
            h = v;
        else
            error(['Provide setRenderFixes with a boolean indicating '...
                   'on/off state for fixes and/or a figure handle']);
        end
        clear v
    end

    ud = get(h, 'UserData');
    n = getNames(h);
    if status
        if isempty(ud)
            set(h, 'UserData',  cell2struct(get(h, n), n, 2));
        end
        % FigureRenderer set to OpenGL to avoid multiscreen issues.
        % DrawMode set to fast to avoid crashes when plotting transparency
        set(h, n{1}, 'OpenGL',...
               n{2}, 'fast');
        warning('mrst:setRenderFixes', ['Default render mode has been set to OpenGL.', ...
            ' This may cause issues with 2D plots such as logarithmic axes']);
    else
        if isempty(ud)
            % In case user clears UserData for whatever reason
            defaults = struct(n{1}, 'painters', ...
                              n{2}, 'normal');
        else
            defaults = ud;
        end
        set(h, n{1}, defaults.(n{1}),...
               n{2}, defaults.(n{2}));
    end
end

function n = getNames(h)
    if h == 0
        n = {'DefaultFigureRenderer', ...
             'DefaultAxesDrawMode'};
    else
        n = {'Renderer', ...
             'DefaultAxesDrawMode'};
    end
end
