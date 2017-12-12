function varargout = plotFractureLines(G,fracture,varargin)
% This function plots the fracture lines using the fracture structure
% returned by processFracture2D.
%
% SYNOPSIS:
%       plotFractureLines(G, fracture)
%       plotFractureLines(G, fracture, 'pn1', pv1, ...)
%   h = plotFractureLines(...)
%
% REQUIRED PARAMETERS:
%
%   G  -        Grid data structure with fractures as assembled by
%               assembleGlobalGrid.
%
%   fracture  - See processFracture2D.
%
% OPTIONAL PARAMETERS:
%
%   show        - Takes string values 'lines' or 'network' to indicate if
%                 each fractures should be coloured using its corresponding
%                 fracture line index or fracture network index,
%                 respectively. Only works if no line numbers are provided
%                 using the second optional argument.
%
%   lineNumbers - Indices of fracture lines to be plotted. If provided, the
%                 above optional parameter is ignored.
%
% RETURNS:
%   h - Handle to polygonal patch structure as defined by function
%       plotFaces.  OPTIONAL.
%
% SEE ALSO:
%   processFracture2D

%{
Copyright 2009-2015: TU Delft and SINTEF ICT, Applied Mathematics.

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

opt = struct('show'       ,  'lines', ...
    'lineNumbers',  []     );
opt = merge_options(opt, varargin{:});

dc = zeros(numel(fracture.lines),3);

%
if isempty(opt.lineNumbers)
    if strcmp(opt.show,'lines')
        dc = rand(numel(fracture.lines),3);
    elseif strcmp(opt.show,'network')
        dn = rand(numel(fracture.network),3);
        for i = 1:numel(fracture.network)
            lines = fracture.network(i).lines;
            dc(lines,:) = repmat(dn(i,:),numel(lines),1);
        end
    end
else
    dc = rand(numel(opt.lineNumbers),3);
end

%

h = plotGrid(G,'FaceColor','w','EdgeAlpha',0.07);
ax = gca;
hold on

%

if ~isempty(opt.lineNumbers)
    lines = opt.lineNumbers;
    for i = 1:numel(lines)
        x = [fracture.lines(lines(i)).endp(1);fracture.lines(lines(i)).endp(3)];
        y = [fracture.lines(lines(i)).endp(2);fracture.lines(lines(i)).endp(4)];
        line(x,y,'Color',dc(i,:),'Parent',ax,'LineWidth',1.5);
        %     text(mean(x),mean(y),num2str(i)); % show line numbering
    end
else
    for i = 1:numel(fracture.lines)
        x = [fracture.lines(i).endp(1);fracture.lines(i).endp(3)];
        y = [fracture.lines(i).endp(2);fracture.lines(i).endp(4)];
        line(x,y,'Color',dc(i,:),'Parent',ax,'LineWidth',1.5);
        %     text(mean(x),mean(y),num2str(i)); % show line numbering
    end
end

if nargout > 0, varargout{1} = h; end

return