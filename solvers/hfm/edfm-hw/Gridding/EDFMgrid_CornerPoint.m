function [G, fracplanes] = EDFMgrid_CornerPoint(G, fracplanes, varargin)
% EDFM for Corner Point grids (Eclipse/Petrel/GRDECL)
% 
% Similar to EDFMgrid but works with Corner Point structure
% Each cell can be an irregular hexahedron with 8 corners
%
% Example:
%   grdecl = readGRDECL('mymodel.grdecl');
%   G = processGRDECL(grdecl);
%   G = computeGeometry(G);
%   [G, frac] = EDFMgrid_CornerPoint(G, fracplanes);

if nargout ~= 2
    error('Output strictly needs to be both G and fracplanes')
end

t1 = clock;

%% Options
opt = struct('Tolerance', 1e-5, ...
             'plotgrid', false, ...
             'fracturelist', 1:length(fracplanes), ...
             'Verbose', true, ...
             'lowesttol', 1e-15);
         
opt = merge_options(opt, varargin{:});
verbose = opt.Verbose;
fracturelist = opt.fracturelist;

%% Rearrange G structure
if ~isfield(G, {'Matrix', 'FracGrid'})
    Gtemp = G;
    Gtemp.Matrix = G;
    Gtemp.FracGrid = struct;
    Gtemp.nnc = [];
    G = Gtemp;    
end

%% Add processed field
if ~isfield(fracplanes, 'processed')
    for i = 1:length(fracplanes)
        fracplanes(i).processed = false;
    end
end

%% Remove preprocessed fractures
templist = [];
for i = fracturelist
    if fracplanes(i).processed
        dispif(verbose, ['Fracplane #', num2str(i), ' has been preprocessed before.']);
        continue;
    end
    templist = [templist, i];
end
fracturelist = templist;

if isempty(fracturelist)
    dispif(verbose, 'All fractures have been preprocessed.');
    return;
end

%% Calculate plane normals
if ~isfield(fracplanes, 'normal')
    for i = 1:numel(fracplanes)
        points = [fracplanes(i).points; fracplanes(i).points(1,:)];
        diffp = diff(points, 1);
        normal = cross(diffp(1,:), diffp(2,:));
        fracplanes(i).normal = normal / norm(normal);
    end
end

%% Extract Corner Point geometry for all cells
Gm = G.Matrix;
nummatcells = Gm.cells.num;

dispif(verbose, 'Extracting Corner Point geometry...\n');

% Extract nodes for all cells
[cn, cpos] = gridCellNodes(Gm, 1:nummatcells);
mcellnodes_par = cell(nummatcells, 1);

for i = 1:nummatcells
    mcellnodeind = cn(cpos(i):(cpos(i+1)-1));
    mcellnodes = Gm.nodes.coords(mcellnodeind, :);
    mcellnodes_par{i} = mcellnodes;
end

% Calculate gaugetol from largest cell
[~, mindex] = max(Gm.cells.volumes);
mcellnodes = mcellnodes_par{mindex};
gaugetol = norm(max(mcellnodes) - min(mcellnodes));

%% Process fractures
Fstorecell = cell(length(fracturelist), 1);
lenfraclist = num2str(length(fracturelist));

for k = 1:length(fracturelist)
    tol = opt.Tolerance;
    i = fracturelist(k);
    
    points = fracplanes(i).points;
    planenormal = fracplanes(i).normal;
    planecenter = sum(points) / size(points, 1);
    aperture = fracplanes(i).aperture;
    
    dispif(verbose, ['Processing ', num2str(k), ' out of ', lenfraclist, ' fractures.\n']);
    
    tightentolerance = true;
    iteration = 0;
    
    while tightentolerance
        iteration = iteration + 1;
        assert(tol > opt.lowesttol, 'Tolerance too tight.');
        
        if iteration > 1
            dispif(verbose, ['  Iteration ', num2str(iteration), ' with tol=', num2str(tol), '\n']);
        end
        
        %% Quick filter
        planepoint = points(1, :);
        temp = bsxfun(@minus, points, planecenter);
        gaugeradius = realsqrt(max(sum(temp.^2, 2)));
        
        mcellsintersect = false(nummatcells, 1);
        
        parfor j = 1:nummatcells
            mcellnodes = mcellnodes_par{j};
            
            % Check intersection with infinite plane
            [xtruth, xtype] = checkplaneAABBintersect_CP(mcellnodes, planepoint, planenormal, tol);
            
            if xtruth && (strcmp(xtype, 'interior') || strcmp(xtype, 'face'))
                % Circular filter
                circtruth = circularfilter(mcellnodes, planenormal, points, tol, ...
                    'gaugetol', gaugetol, 'gaugeradius', gaugeradius, 'center', planecenter);
                if circtruth
                    mcellsintersect(j) = true;
                end
            end
        end
        
        possiblemcells = find(mcellsintersect);
        dispif(verbose, ['  ', num2str(length(possiblemcells)), '/', num2str(nummatcells), ' cells selected.\n']);
        
        if isempty(possiblemcells)
            warning('No cells intersect fracture %d', i);
            break;
        end
        
        %% Rigorous intersection detection - Corner Point Specific
        Nm = length(possiblemcells);
        
        mcells = -1 * ones(Nm, 1);
        area = -1 * ones(Nm, 1);
        typenum = -1 * ones(Nm, 1);
        fraccellpoints = cell(Nm, 1);
        centroidside = false(Nm, 1);
        
        % Prepare data for parfor
        mcellnodes_par_i = cell(Nm, 1);
        for j = 1:Nm
            mcellnodeind = cn(cpos(possiblemcells(j)):(cpos(possiblemcells(j)+1)-1));
            mcellnodes = Gm.nodes.coords(mcellnodeind, :);
            mcellnodes_par_i(j) = {mcellnodes};
        end
        
        parfor j = 1:Nm
            mcellnodes = mcellnodes_par_i{j};
            
            % Intersect fracture plane with Corner Point hexahedron
            [xtruth, xarea, xtype, xlocation, fracpts] = pebiHexahedronIntersect_CP(...
                points, mcellnodes, tol);
            
            if xtruth && strcmp(xtype, 'polygon')
                mcells(j) = possiblemcells(j);
                area(j) = xarea;
                
                switch xlocation
                    case 'none'
                        xlocation = 0;
                    case 'boundary'
                        xlocation = 1;
                        centroid = sum(mcellnodes) / size(mcellnodes, 1);
                        dir = dot(centroid - points(1,:), planenormal);
                        dir = dir / abs(dir);
                        centroidside(j) = dir > 0;
                    case 'interior'
                        xlocation = 2;
                        centroidside(j) = true;
                end
                
                typenum(j) = xlocation;
                fraccellpoints{j} = fracpts;
            end
        end
        
        %% Filter and sort
        intersected = ~cellfun('isempty', fraccellpoints);
        mcells = [mcells(intersected & centroidside); mcells(intersected & ~centroidside)];
        area = [area(intersected & centroidside); area(intersected & ~centroidside)];
        typenum = [typenum(intersected & centroidside); typenum(intersected & ~centroidside)];
        
        fraccellpoints = fraccellpoints(intersected & centroidside);
        
        if isempty(fraccellpoints)
            warning('No valid intersection for fracture %d', i);
            break;
        end
        
        %% Construct V and C
        V = vertcat(fraccellpoints{:});
        
        C = cellfun(@(c) 1:size(c,1), fraccellpoints, 'UniformOutput', false);
        for j = 2:size(C, 1)
            addTo = C{j-1}(end);
            C{j} = C{j} + addTo;
        end
        
        %% Create fracture grid
        try
            Fgrid = fractureplanegeneralgrid(V, C, points, planenormal, aperture, tol);
            
            % Matrix-fracture connection
            Fgrid.matrix_connection.cells = mcells;
            Fgrid.matrix_connection.area = area;
            Fgrid.matrix_connection.type = typenum;
            
            % Check volume
            targetvol = convexpolygonarea(points, opt.Tolerance) * aperture;
            totcellvol = sum(Fgrid.cells.volumes);
            relerror = abs(targetvol - totcellvol);
            
            dispif(verbose, ['  Target vol: ', num2str(targetvol), ...
                           ', Computed: ', num2str(totcellvol), ...
                           ', Error: ', num2str(relerror), '\n']);
            
            tightentolerance = ~(relerror < opt.Tolerance);
            
            if ~tightentolerance
                Fstorecell{k} = Fgrid;
            end
            
       catch ME
    warning(ME.identifier, '%s', ME.message);
end

        
        tol = tol / 10;
        
            if  iteration > 10
                warning('Max iterations exceeded!');
                break;
            end
    end
end

%% Append to FracGrid
FracGrid_new = G.FracGrid;
cstart = G.cells.num + 1;
fstart = G.faces.num + 1;
nstart = G.nodes.num + 1;
lastfrac = length(fieldnames(G.FracGrid));

for j = 1:length(fracturelist)
    if isempty(Fstorecell{j})
        continue;
    end
    
    i = fracturelist(j);
    fieldname = ['Frac', num2str(j + lastfrac)];
    
    FracGrid_new.(fieldname) = Fstorecell{j};
    FracGrid_new.(fieldname).cells.start = cstart;
    FracGrid_new.(fieldname).faces.start = fstart;
    FracGrid_new.(fieldname).nodes.start = nstart;
    
    cstart = cstart + FracGrid_new.(fieldname).cells.num;
    fstart = fstart + FracGrid_new.(fieldname).faces.num;
    nstart = nstart + FracGrid_new.(fieldname).nodes.num;
    
    FracGrid_new.(fieldname).rock.perm = ...
        ones(FracGrid_new.(fieldname).cells.num, 1) * fracplanes(i).perm;
    FracGrid_new.(fieldname).rock.poro = ...
        ones(FracGrid_new.(fieldname).cells.num, 1) * fracplanes(i).poro;
    
    FracGrid_new.(fieldname).points = fracplanes(i).points;
    FracGrid_new.(fieldname).fracplanenumber = i;
    FracGrid_new.(fieldname).matrixnnc = false;
    FracGrid_new.(fieldname).fracgridnnc = [];
    FracGrid_new.(fieldname).gridtype = FracGrid_new.(fieldname).type;
    FracGrid_new.(fieldname) = rmfield(FracGrid_new.(fieldname), 'type');
    
    fracplanes(i).fracgrid = j + lastfrac;
    fracplanes(i).processed = true;
end

G.FracGrid = FracGrid_new;

%% Plot grid
if opt.plotgrid
    figure;
    plotGrid(G.Matrix, 'facealpha', 0);
    for i = 1:numel(fieldnames(G.FracGrid))
        fieldname = ['Frac', num2str(i)];
        if ismember(G.FracGrid.(fieldname).fracplanenumber, opt.fracturelist)
            facecolor = 'y';
        else
            facecolor = 'r';
        end
        plotGrid(G.FracGrid.(fieldname), 'facealpha', 0.75, 'Facecolor', facecolor);
    end
    view(15, 20);
    axis equal tight;
end

%% Regenerate global grid
Gtemp = G.Matrix;
Gtemp.FracGrid = G.FracGrid;
Gtemp.nnc = G.nnc;
G = assembleGlobalGrid(Gtemp);

t2 = clock;
e = etime(t2, t1);
dispif(verbose, ['EDFM grid assembly took ', num2str(e), ' seconds.\n']);

end

%% ========================================================================
%% Helper Functions - Corner Point Specific
%% ========================================================================

function [xtruth, xtype] = checkplaneAABBintersect_CP(hexnodes, planepoint, planenormal, tol)
% Check if Corner Point hexahedron intersects infinite plane
    
    xtruth = false;
    xtype = 'none';
    
    % Distance of all points to plane
    dists = (hexnodes - planepoint) * planenormal';
    
    % If all points on one side, no intersection
    if all(dists > tol) || all(dists < -tol)
        return;
    end
    
    xtruth = true;
    
    % Determine intersection type
    if any(abs(dists) < tol)
        xtype = 'face';  % plane passes through face
    else
        xtype = 'interior';  % plane passes through interior
    end
end

function [xtruth, xarea, xtype, xlocation, clipped] = pebiHexahedronIntersect_CP(...
    fracPolygon, hexnodes, tol)
% Rigorous intersection of fracture plane with Corner Point hexahedron
% Uses Sutherland-Hodgman algorithm
    
    xtruth = false;
    xarea = 0;
    xtype = 'none';
    xlocation = 'none';
    clipped = [];
    
    % Clip polygon with 6 faces of hexahedron
    clipped = clipPolygonByHexahedron(fracPolygon, hexnodes, tol);
    
    if isempty(clipped) || size(clipped, 1) < 3
        return;
    end
    
    % Calculate area
    xarea = convexpolygonarea(clipped, tol);
    
    if xarea < tol^2
        return;
    end
    
    xtruth = true;
    xtype = 'polygon';
    
    % Determine location
    if size(clipped, 1) == size(fracPolygon, 1)
        xlocation = 'interior';
    else
        xlocation = 'boundary';
    end
end

function clipped = clipPolygonByHexahedron(polygon, hexnodes, tol)
% Clip polygon with 6 faces of hexahedron
% Corner Point node order typically:
% Bottom: 1-2-3-4, Top: 5-6-7-8
    
    clipped = polygon;
    
    nNodes = size(hexnodes, 1);
    
    % Define 6 faces of hexahedron
    if nNodes == 8
        % Standard hexahedron
        faces = {
            [1, 4, 3, 2],  % bottom
            [5, 6, 7, 8],  % top
            [1, 2, 6, 5],  % front
            [4, 8, 7, 3],  % back
            [1, 5, 8, 4],  % left
            [2, 3, 7, 6]   % right
        };
    else
        % Non-standard number of nodes
        warning('Non-standard hexahedron with %d nodes', nNodes);
        return;
    end
    
    for f = 1:6
        if size(clipped, 1) < 3
            break;
        end
        
        faceIdx = faces{f};
        
        % Check valid indices
        if any(faceIdx > nNodes)
            continue;
        end
        
        faceNodes = hexnodes(faceIdx, :);
        
        % Calculate face normal
        v1 = faceNodes(2,:) - faceNodes(1,:);
        v2 = faceNodes(3,:) - faceNodes(1,:);
        normal = cross(v1, v2);
        
        if norm(normal) < tol
            continue;
        end
        
        normal = normal / norm(normal);
        
        % Clip polygon with this plane
        clipped = clipPolygonByPlane(clipped, faceNodes(1,:), normal, tol);
    end
end

function clipped = clipPolygonByPlane(polygon, planepoint, planenormal, tol)
% Sutherland-Hodgman algorithm
    
    if isempty(polygon) || size(polygon, 1) < 3
        clipped = [];
        return;
    end
    
    output = [];
    n = size(polygon, 1);
    
    for i = 1:n
        curr = polygon(i, :);
        next = polygon(mod(i, n) + 1, :);
        
        dcurr = dot(curr - planepoint, planenormal);
        dnext = dot(next - planepoint, planenormal);
        
        if dcurr >= -tol
            if dnext >= -tol
                output = [output; next];
            else
                t = dcurr / (dcurr - dnext + eps);
                isect = curr + t * (next - curr);
                output = [output; isect];
            end
        else
            if dnext >= -tol
                t = dcurr / (dcurr - dnext + eps);
                isect = curr + t * (next - curr);
                output = [output; isect; next];
            end
        end
    end
    
    if ~isempty(output)
        output = uniquetol(output, tol, 'ByRows', true);
    end
    
    clipped = output;
end