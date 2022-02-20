function varargout = mrstFigure(varargin)
%Create or select figure with convenient default options for 3d plotting
%
% SYNOPSIS:
%   h = mrstFigure();
%   mrstFigure(h);
%   mrstFigure(h, 'background', 'fast');
%
% DESCRIPTION:
%   This function is intended to act as figure() when the figure is
%   intended to plot 3D datasets. This allows for easier syntax and ensures
%   that reasonable default parameters can be set. For instance, the
%   default renderer will be set to OpenGL, which will avoid certain
%   crashes when plotting programatically.
%
% OPTIONAL PARAMETERS:
%   Supply as keywords, i.e. mrstFigure('keyword').
%
%   mrstFigure(h, 'background') - set current figure to h, without making
%   the window pop up.
%
%   mrstFigure(h, 'fast') - set DrawMode to fast for faster drawing, but
%   somewhat uglier appearance in presence of transparency. Note that this
%   can prevent some crashes.
%
%   mrstFigure(h, 'matlabrotate') - Use standard MATLAB rotate function
%   which may be slower.
%
%   mrstFigure(h, 'painters') - Use default renderer instead of OpenGL.
%
% RETURNS:
%   h - handle to current figure.
%
% EXAMPLE:
%    G = cartGrid([10 10 1]);
%    mrstFigure();
%    plotGrid(G);
%
% SEE ALSO:
%   `figure`, `fastRotateButton`, `setRenderFixes`

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

    % Check if first input argument is a handle to a valid figure,
    if numel(varargin)
        h = varargin{1};
        varargin = varargin(2:end);
    else
        h = figure();
    end
    
    % Silently set figure
    if any(strcmpi(varargin, 'background')) && validateFunctionHandle(h)
        set(0, 'CurrentFigure', h);
    else
        figure(h);
    end
    
    % Set properties and fast rotation
    if any(strcmpi(varargin, 'fast'))
        set(h, 'DefaultAxesDrawMode', 'fast')
    end
    
    if ~any(strcmpi(varargin, 'painters'))
        set(h, 'Renderer', 'OpenGL')
    end
    
    if ~any(strcmpi(varargin, 'matlabrotate'))
        fastRotateButton(h);
    end
    
    if nargout
        varargout{1} = h;
    end
end

function isFigureHandle = validateFunctionHandle(h)
    if numel(double(h)) == 1 && ishandle(h)...
                             && strcmpi(get(h, 'Type'), 'figure')
        isFigureHandle = true;
    else
        isFigureHandle = false;
    end
end