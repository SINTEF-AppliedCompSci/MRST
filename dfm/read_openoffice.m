function [vertices, edges, box] = read_openoffice(filename, ext, opts)
% Read fractures drawn in Libreoffice, and represent them as points and
% edges.
%
% Copyright (C) 2006-2007 Uni Research AS
% This file is licensed under the GNU General Public License v3.0.

if exist(filename, 'file') ~= 2,
   error(['Fracture input file ''%s'' ', ...
          'does not exist on local system'], filename);
end

% directory of this function contains the stylesheets
dir = fileparts (which (mfilename ()));

% start out with empty sets; all information will be read from the file
vertices = [];

% extract the style file and use it to find out the maximum and minimum
% coordinates of the figure
[stat, res] = ...
   system (sprintf (['unzip -p "%s" styles.xml |', ...
                     'xsltproc --novalid "%s" - > %sbox.mat'], ...
                     filename, fullfile (dir, ['box-', ext, '.xsl']), ...
           [dir, filesep]));

if stat ~= 0,
   error('Failed to convert fracture input file ''%s'':\n\t%s', ...
         filename, res);
end

% create a bounding box from the min./max. coordinates

maxcoords = load (fullfile(dir, 'box.mat'), '-ASCII');
box = [0, 0; maxcoords(1), maxcoords(2)];

% extract the content file and run it through the style sheet to get to
% we use a pipe instead of the xslt() function so we don't have to store
% the entire content file ourself.
system (sprintf (['unzip -p "%s" content.xml | xsltproc --novalid ' ...
    '"%s" - > lines.mat'], filename, fullfile (dir, ['lines-' ext '.xsl'])));

% read the coordinates that was extracted; from these coordinates we can
% create a set of vertices and lines
linecoords = load ('lines.mat', '-ASCII');

% extract the content file and run it through the style sheet to get to
% we use a pipe instead of the xslt() function so we don't have to store
% the entire content file ourselves.
system (sprintf (['unzip -p "%s" content.xml | xsltproc --novalid ' ...
    '"%s" - > type.mat'], filename, fullfile (dir, ['type-' ext '.xsl'])));

% read the type of the lines
type = load ('type.mat', '-ASCII');

% READS THE COLOR OF THE LINES may be useful
%    % extract the content file and run it through the style sheet to get to
%   % we use a pipe instead of the xslt() function so we don't have to store
%   % the entire content file ourself.
system (sprintf (['unzip -p "%s" content.xml | xsltproc --novalid ' ...
    '"%s" - > color.mat'], filename, fullfile (dir, ['color-' ext '.xsl'])));

if isfield(opts,'color')
    c = load ('color.mat', '-ASCII');
    color2tagmap = zeros(numel(c),1);
    for i = 1:numel(c)
        c_str = num2str(c(i));
        l_str = numel(c_str);
        while l_str<6
            c_str = ['0',c_str];
            l_str = numel(c_str);
        end
        color = hex2dec({c_str(1:2),c_str(3:4),c_str(5:6)})'/128;
        color2tagmap(i) = find(all(bsxfun(@eq,color,opts.color),2));
    end
    t = color2tagmap(type);

else
    % currently, all fractures are of type 0. we could map each color to its own
    % type and use that information to decide which multiplier it should be
    % assigned.
    t = type;
end



% create vertices and edges from the start and end coordinates
% respectively. each row in the file is one line with four values on it
edges = zeros (size (linecoords, 1), 3);
for i = 1:size (linecoords, 1)
    [vertices, a] = ...
        add_point (vertices, [linecoords(i, 1), linecoords(i, 2)], box, opts);
    [vertices, b] = ...
        add_point (vertices, [linecoords(i, 3), linecoords(i, 4)], box, opts);
    edges(i, :) = [a, b, t(i)];
end;

% remove temporary file after it is loaded
system ('rm lines.mat');
system ('rm box.mat');
system ('rm type.mat');
system ('rm color.mat');

% second coordinate is swapped when reading from OpenOffice; we need to
% turn it upside-down.
y = 2; min = 1; max = 2;
vertices(:, y) = box(max, y) - (vertices(:, y) - box(min, y)) + box(min, y);
