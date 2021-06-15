function printStateFunctionGroupingTikz(g, varargin)
%Undocumented Utility Function

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

    if mod(numel(varargin), 2) == 1
        p = varargin{1};
        varargin = varargin(2:end);
    else
        p = nan;
    end
    opt = struct('edgeDraw',    true, ...
                 'drawGroups',  true, ...
                 'inner',       'tree layout, grow = right', ...
                 'outer',       'tree layout', ...
                 'file',        []);
    opt = merge_options(opt, varargin{:});
    if isempty(opt.file)
        fp = 1;
    else
        fp = fopen(opt.file,'w');
    end
    % Tikz version of plotStateFunctionGrouping
    outer = opt.outer;
    inner = opt.inner;
    node_dist = 3;
    nodes = g.Nodes{:, 1};
    edges = g.Edges{:, 1};
    n = numel(nodes);
    nodeNames = cell(1, n);
    groupNames = cell(1, n);
    hasHandle = ishandle(p);
    if hasHandle
        nodeNames = get(p, 'NodeLabel');
        for i = 1:numel(nodeNames)
            nn = nodeNames{i};
            if checkLaTeX(nn)
                nodeNames{i} = sprintf('$%s$', nn);
            end
        end
    end
    for i = 1:n
        node = nodes{i};
        sep = regexp(node,'\.','split');
        if numel(sep) == 1
            nodeName = sep{1};
            groupName = 'state';
        else
            nodeName = sep{2};
            groupName = sep{1};
        end
        if ~hasHandle
            nodeNames{i} = nodeName;
        end
        groupNames{i} = groupName;
    end
    uniqueGroups = uniqueStable(groupNames);
    isState = strcmp(uniqueGroups, 'state');
    if any(isState)
        tmp = uniqueGroups{1};
        uniqueGroups{1} = 'state';
        uniqueGroups{isState} = tmp;
    end
    ng = numel(uniqueGroups);
    
    fprintf(fp, ['\\documentclass[tikz,border=5pt]{standalone}\n', ...
            '\\usepackage{tikz}\n', ...
            '\\usetikzlibrary{graphdrawing.layered}\n', ...
            '\\usetikzlibrary{positioning,shapes}\n', ...
            '\\usetikzlibrary{arrows}\n', ...
            '\\usetikzlibrary{graphs}\n', ...
            '\\usetikzlibrary{graphdrawing}\n', ...
            '\\usetikzlibrary{backgrounds}\n', ...
            '\\usetikzlibrary{svg.path}\n', ...
            '\\usetikzlibrary{bending}\n', ...
            '\\usetikzlibrary{shapes.geometric,arrows.meta,decorations.markings}\n', ...
            '\\usegdlibrary{layered, trees, force, circular, phylogenetics}\n', ...
            '\\begin{document}\n']);
    
    fprintf(fp, '\\tikzstyle{propbox}=[rounded rectangle, draw = black]\n');
    fprintf(fp, '\\tikzstyle{propedge}=[>={Stealth[round,sep,bend]}, line width = 1pt, draw = black!30, line width = 1pt, draw = black!30, -Stealth]\n');
    fprintf(fp, '\\tikzstyle{groupbox}=[font=\\bfseries \\Large, rounded corners, opacity=0.4]\n');
    
    colors = [228, 26, 28; ...
              55, 126, 184;...
              77, 175, 74; ...
              152, 78, 163; ...
              255, 127, 0]./255; % http://colorbrewer2.org/#type=qualitative&scheme=Set1&n=5
    colors = repmat(colors, ng, 1);
    colors = colors(1:ng, :);
    for i = 1:ng
        c1 = colors(i, :);
        c2 = brighten(c1, 0.5);
        fprintf(fp, '\\tikzstyle{%s}=[propbox, fill={%s}]\n', uniqueGroups{i}, formatRGB(c2));
        fprintf(fp, '\\tikzstyle{%sEdge}=[propedge, draw={%s}]\n', uniqueGroups{i}, formatRGB(c1));
    end
    fprintf(fp, '\\tikz[]{\n');
    fprintf(fp, '\\graph[%s, nodes={propbox}, edges={propedge}, node distance = %dcm]{\n', outer, node_dist);
    active = false(size(edges, 1), 1);
    doSub = opt.drawGroups;
    
    for gno = 1:ng
        group = uniqueGroups{gno};
        if doSub
            fprintf(fp, '%s[font=\\bfseries \\Large] // [%s,  edges={%sEdge}] {\n', group, inner, group);
        end
        sub = find(strcmp(groupNames, group));
        % Draw nodes
        for index = 1:numel(sub)
            i = sub(index);
            fprintf(fp, '%d [%s, as=%s], \n', i, group, nodeNames{i});
        end
        act = drawEdges(fp, nodes, edges, groupNames, group, [], []);
        active = active | act;
        if doSub
            fprintf(fp, '},\n');
        end
    end
    % Draw edges
    edgecolors = colors;
    if ~opt.edgeDraw
        drawEdges(fp, nodes, edges(~active, :), groupNames, [], 'densely dashed, opacity = 0.5', edgecolors, uniqueGroups);
    end
    fprintf(fp, '};\n');
    fprintf(fp, '\\begin{scope}[on background layer]\n');
    if opt.edgeDraw
        drawEdges(fp, nodes, edges(~active, :), groupNames, [], 'propedge, densely dashed, opacity = 0.5', edgecolors, uniqueGroups, true);
    end
    % Draw background boxes
    if doSub
        for gno = 1:ng
            group = uniqueGroups{gno};
            fprintf(fp, '\t\\draw[%s, groupbox]\n', group);
            fprintf(fp, '\t(%s.north east) rectangle (%s.south west);\n', group, group);
        end
    end
    fprintf(fp, '\\end{scope}\n');
    fprintf(fp, '}\n\\end{document}\n');
    if ~isempty(opt.file)
        fclose(fp);
    end
end

function istex = checkLaTeX(x)
    if numel(x) <= 3
        istex = true;
    elseif any(x == '{' | x == '}' | x == '_' | x == '\')
        istex = true;
    else
        istex = false;
    end
end

function drawn = drawEdges(fp, nodes, edges, groups, group, stylearg, colors, allgroups, asSep)
    if nargin < 9
        asSep = false;
    end
    ne = size(edges, 1);
    drawn = false(ne, 1);
    for i = 1:ne
        edge = edges(i, :);
        start = find(strcmp(nodes, edge{1}));
        stop = find(strcmp(nodes, edge{2}));
        if ~isempty(group)
            if ~(strcmp(groups{stop}, group) && strcmp(groups{start}, group))
                continue
            end
        end
        drawn(i) = true;
        
        arg = stylearg;
        if ~isempty(colors)
            c = colors(strcmp(groups{start}, allgroups), :);
            arg = [arg, ',', sprintf('draw={%s}', formatRGB(c))];
        end
        if isempty(arg)
            arg = '';
        else
            arg = ['[', arg, ']'];
        end
        if asSep
            fprintf(fp, '\\draw%s (%d) -> (%d);\n', arg, start, stop);
        else
            fprintf(fp, '%d -> %s %d, \n', start, arg, stop);
        end
    end
end

function s = formatRGB(rgb)
    rgb = ceil(rgb*255);
    s = sprintf('rgb,255:red,%d; green,%d; blue,%d', rgb);
end
