function varargout = plotWell(G, W, varargin)
%Plot well trajectories into current axes.
%
% SYNOPSIS:
%   [htop, htext, hs, hline] = plotWell(G, W)
%   [htop, htext, hs, hline] = plotWell(G, W, 'pn1', pv1, ...)
%
% PARAMETERS:
%   G       - Grid data structure
%   W       - Well data structure
%
% KEYWORD ARGUMENTS:
%  'radius' - Scaling factor for the width of the well path. 
%             Default value: 1.0
%
%  'height' - Height above top reservoir contact at which the well should
%             stop and symbol be drawn. 
%             Default value: height = 5
%
%  'color'  - Colour with which the well path should be drawn.  Possible
%             values described in `plot`.
%             Default value: color = 'r'
%
%  'color2' - Second color. Will be used for producer wells. If not
%             specified, all wells will have the same color and the `color`
%             argument will be used.
%
%  'cylpts' - Number of segments to use about the well bore. A higher value
%             produces more smoothly looking well trajectories at the
%             expense of more costly plotting.
%             Default value: cylpts = 10
%
%  'fontsize' - The size of the font used for the label texts.
%               Default value: fontsize = 16
%
%  'ambstr' - The ambient strength of the well cylinder.
%             Default value: ambstr = 0.8
%
%  'linewidth' - The width of the line above the well.
%                Default value: linewidth = 2
%
% RETURNS:
%   htop  - Graphics handles to all well heads and heels.
%   htext - Graphics handles to all rendered well names.
%   hs    - Graphics handles to all well paths.
%   hline - Graphics handles to all lines between well head and names
%
% EXAMPLE:
%   G = computeGeometry(cartGrid([10, 1, 1]));
%   rock = makeRock(G, 1, 1);
%   W = addWell([], G, rock, 1)
%   plotWell(G, W, 'color', 'k');
%
% NOTE:
%   Currently, this function only supports three-dimensional grids.
%
% SEE ALSO:
%   `addWell`, `delete`, `patch`, `incompTutorialWells`.

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


assert(G.griddim==3,'plotWell not implemented for 2D grids');

prm  = struct('radius',    1,  ...
              'height',    5,  ...
              'color',     [],  ...
              'color2',    [], ...
              'cylpts',    10,  ...
              'fontsize',  16,  ...
              'ambstr',    0.8,  ...
              'linewidth', 2);
prm = merge_options(prm, varargin{:});

dohold = ishold;
hold on;

nW    = numel(W);
% Number of digits required to textually represent all wells on the form
% W%d (e.g. 'W03') with the same number of digits for all wells.
%
nWd = fix(log10(nW)) + 1;

htop  = zeros([1, 2*nW]);
htext = zeros([1,   nW]);
hs    = zeros([1,   nW]);
hline = zeros([1,   nW]);

% Common face/vertex connection table for plotting well paths by means of
% PATCH.
faces = @(sz,i,j) [sub2ind(sz, i  , j  ), ... % lower left  vertex
                   sub2ind(sz, i+1, j  ), ... % lower right vertex
                   sub2ind(sz, i+1, j+1), ... % upper right vertex
                   sub2ind(sz, i  , j+1)];    % upper left  vertex
if isCoarseGrid(G)
    zcoord = G.parent.nodes.coords(:, 3);
else
    zcoord = G.nodes.coords(:, 3);
end
ztop =  min(zcoord) - prm.height;
for w = 1 : nW
   if W(w).sign < 0 && ~isempty(prm.color2)
       color = prm.color2;
   elseif ~isempty(prm.color)
       color = prm.color;
   else
       color = 'r';
   end
   % Plot cylinder
   c = G.cells.centroids(W(w).cells, :);

   dir = W(w).dir;
   if numel(W(w).cells) == 1
      dir = 's';
   end

   % Generate a cylinder rotated according to major direction
   [XW,YW,ZW] = deal([]);
   for i=1:numel(dir)
      switch lower(dir(i))
         case 's'
            [xw, yw, zw] = sphere(prm.cylpts);
            cno = W(w).cells(i);
            fno = G.cells.facePos(cno):G.cells.facePos(cno+1)-1;
            area = max(G.faces.areas(G.cells.faces(fno)));
            rl = prm.radius*0.25*sqrt(area);
            rv = prm.radius*0.25*G.cells.volumes(cno)/area;
            xw = rl*xw + repmat(c(i,1), (prm.cylpts + 1)*[1 1]); XW = [XW; xw]; %#ok
            yw = rl*yw + repmat(c(i,2), (prm.cylpts + 1)*[1 1]); YW = [YW; yw]; %#ok
            zw = rv*zw + repmat(c(i,3), (prm.cylpts + 1)*[1 1]); ZW = [ZW; zw]; %#ok
         otherwise
            continue
      end
   end

   for i=1:numel(dir)-1
      cno = W(w).cells(i);
      fno = G.cells.facePos(cno):G.cells.facePos(cno+1)-1;
      area = max(G.faces.areas(G.cells.faces(fno)));
      rl = prm.radius*0.25*sqrt(area);
      rv = prm.radius*0.25*G.cells.volumes(cno)/area;
      switch lower(dir(i))
         case 'x'
            R = [rv; rv];
            [zw, yw, xw] = cylinder(R, prm.cylpts);
            xw = 0*xw;
         case 'y'
            R = [rv; rv];
            [xw, zw, yw] = cylinder(R, prm.cylpts);
            yw = 0*yw;
         case 'z'
            R = [rl; rl];
            [xw, yw, zw] = cylinder(R, prm.cylpts);
            zw = 0*zw;
         case 's'
            continue;
         otherwise
            R = [rl; rl];
            [xw, yw, zw] = cylinder(R, prm.cylpts);
            zw = rv*(zw-0.5);
      end
      xw = xw + repmat(c(i:i+1,1), [1, prm.cylpts + 1]); XW = [XW; xw]; %#ok
      yw = yw + repmat(c(i:i+1,2), [1, prm.cylpts + 1]); YW = [YW; yw]; %#ok
      zw = zw + repmat(c(i:i+1,3), [1, prm.cylpts + 1]); ZW = [ZW; zw]; %#ok
   end
   %zw(1,:) = zw(1,:) - prm.height;

   % Draw top
   nodes = [XW(1,:); YW(1,:); ZW(1,:)] .';
   htop(2*(w-1) + 1) = patch('Faces',  (1 : numel(nodes)/3), ...
                             'Vertices',  nodes,             ...
                             'EdgeColor', 'k',               ...
                             'FaceColor', color,         ...
                             'AmbientStrength', prm.ambstr);

   % Draw bottom
   nodes = [XW(end,:); YW(end,:); ZW(end,:)] .';
   htop(2*w) = patch('Faces',  (1 : numel(nodes)/3), ...
                     'Vertices', nodes,              ...
                     'EdgeColor', 'k',               ...
                     'FaceColor', color,         ...
                     'AmbientStrength', prm.ambstr);

   if isfield(W(w), 'name') && ischar(W(w).name)
      name = escape_for_LaTeX_interpreter(W(w).name);
   else
      name = sprintf('W$_{%0*d}$', nWd, w);
   end
   if prm.fontsize > 0
       htext(w) = text(c(1,1),                   ...
                       c(1,2),                   ...
                       ztop,                     ...
                       name, 'FontSize', prm.fontsize,  ...
                       'VerticalAlignment', 'bottom', ...
                       'interp', 'latex', 'Color', color);
   end
   sz    = size(XW);
   [I,J] = ndgrid(1 : sz(1)-1, 1 : sz(2)-1);
   hs(w) = patch('Faces'          , faces(sz, I(:), J(:)), ...
                 'Vertices'       , [XW(:), YW(:), ZW(:)], ...
                 'AmbientStrength', prm.ambstr,            ...
                 'FaceColor'      , color,             ...
                 'EdgeColor'      , 'none');

   % Plot a line to somewhere above the model
   c = c([1 1:end],:);
   c(1,3) = min(zcoord) - prm.height;
   hline(w) = plot3(c(:,1), c(:,2), c(:,3), ...
             'color', color, 'LineWidth', prm.linewidth);

  if nargout > 0
      varargout{1} = htop;
  end

  if nargout > 1
      varargout{2} = htext;
  end

  if nargout > 2
      varargout{3} = hs;
  end

  if nargout > 3
      varargout{4} = hline;
  end
end

if dohold
   hold on;
else
   hold off;
end

%--------------------------------------------------------------------------

function s = escape_for_LaTeX_interpreter(s)
assert (ischar(s), 'Internal Logic Error');

matches = @(re) ~ isempty(regexp(s, re, 'once'));

if matches('_') && ~matches('$')
   % String 's' contains at least one underscore but does *not* contain any
   % '$' characters (which could mean that the name has already been
   % escaped for the LaTeX interpreter of TEXT.)  Escape all '_' characters
   % (i.e., replace by '\_') to satisfy syntax requirements of interpreter.

   s = regexprep(s, '_', '\\_');
end
