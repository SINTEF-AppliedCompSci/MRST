function [pSurfs, pSurfXY, wellbores] = generateHWGridNodes(GV, packed, well, radPara)
% Generate 3D points of all radial HW grid surfaces and 2D planar points
%
% SYNOPSIS:
%   [pSurfs, pSurfXY, wellbores] = generateHWGridNodes(GV, packed, well, radPara)
%
% PARAMETERS:
%  GV           - Unstructured VOI grid structure
%  packed       - Structure of HW information, see 'allInfoOfRegion'
%  well         - Structure of well information, see 
%                 'nearWellBoreModelingGrids'
%  radPara      - Parameters for generating the radial grid, see 
%                 'ReConstructToRadialGrid'
%
% RETURNS:
%  pSurfs     - 3D points of all HW grid surfaces
%  pSurfXY    - pSurfs in xy plane, used to compute the radial
%               transmissibility factors
%  wellbores  - Structure of the casing and screen, used to generate nodes
%               segments for the multi-segment well
%
% SEE ALSO:
%  `HorWellRegion` `generateVOIGridNodes`

    switch radPara.gridType
        case 'pureCircular'
            generator = @(th)CircularRadialPoints(GV, packed, well, radPara, th);
        case 'gradual'
            generator = @(th)GradualRadialPoints(GV, packed, well, radPara, th);
        otherwise
            error([mfilename, ': Unknown radial grid type'])
    end
    
    ns = well.segmentNum;
    [pSurfs, pSurfXY, wellbores] = ...
        arrayfun(generator, (1:ns+1)', 'UniformOutput', false);
end    

% -------------------------------------------------------------------------
function [p, pxy, wellbore] = CircularRadialPoints(GV, packed, well, radPara, ii)
% Generate pure circular radial points for single well node ii
    % Assign parameters
    rM     = radPara.maxRadius;
    nRT    = radPara.nRadCells;
    bdn    = packed.bdyNodes{ii};
    vxID   = packed.vertexID;
    origin = well.trajectory(ii, :);
    rW     = well.radius(ii);
    
    % Coordinate transformation for outer boundary nodes and well center
    [pob0, origin, T, R] = convertToXYPlane(GV.nodes.coords(bdn,:), vxID, origin);
    pob = bsxfun(@minus, pob0, origin);
    pob = pob(:,1:2);
    
    % Theta and radii of the boundary nodes
    [theta, rob] = cart2pol(pob(:,1), pob(:,2));
    assert(all(rob > rM), ['Max radius of the radial grid is greater '...
        ' than the boundary radius, please reduce the value']);
    
    % Generate circular radial points
    rr = exp( linspace(log(rW), log(rM), nRT(1)) );
    [RR, THETA] = meshgrid(rr, theta);
    [px, py] = pol2cart(THETA(:), RR(:));
    
    % Add the boundary points
    px = [px; pob(:,1)];
    py = [py; pob(:,2)];
    
    % Return the planner points
    pxy = [px(:), py(:)];
    
    % Map back to original coordinate
    px = px(:) + origin(1);
    py = py(:) + origin(2);
    
    % Get pz by interpolation
    xx = [pob0(:,1); origin(1)];
    yy = [pob0(:,2); origin(2)];
    zz = [pob0(:,3); origin(3)];
    pz = griddata(xx, yy, zz, px, py);
    pz(1:length(bdn)) = pob0(:,3);
    
    % Return the 3D points
    p = [px, py, pz];
    p = convertTo3DPlane(p, T, R);
    
    % Store borewall and screen nodes, prepared for multi-segment well
    nA = length(bdn);
    pW = p(1:nA, :);
    wellbore.wall.radius = rW;
    wellbore.wall.coords = pW;
    wellbore.reservoirCells = (1:nA)';
    if isfield(well, 'screenRadius')
        rS  = well.screenRadius(ii);
        assert(rS < rW, ...
            'The screen radius must be less than the casing radius!')
        [pxS, pyS] = pol2cart(theta, rS);
        pxS = pxS + origin(1);
        pyS = pyS + origin(2);
        pzS = griddata(xx, yy, zz, pxS, pyS);
        pS = [pxS, pyS, pzS];
        pS = convertTo3DPlane(pS, T, R);
        wellbore.screen.radius = rS;
        wellbore.screen.coords = pS;
    end
end

% -------------------------------------------------------------------------
function [p, pxy, wellbore] = GradualRadialPoints(GV, packed, well, radPara, ii)
% Generate gradual radial points for single well node ii
    % Assign parameters
    boxRatio  = radPara.boxRatio;
    nRT       = radPara.nRadCells;
    pDMult    = radPara.pDMult;
    offCenter = radPara.offCenter;
    bdn       = packed.bdyNodes{ii};
    vxID      = packed.vertexID;
    origin    = well.trajectory(ii, :);
    rW        = well.radius(ii);
    
    % Coordinate transformation for outer boundary nodes and well center
    [pob0, origin, T, R] = convertToXYPlane(GV.nodes.coords(bdn,:), vxID, origin);
    pob = bsxfun(@minus, pob0, origin);
    pob = pob(:,1:2);
    
    % Rectangular box size
    a = boxRatio(1) * (max(pob(:,1)) - min(pob(:,1)));
    b = boxRatio(2) * (max(pob(:,2)) - min(pob(:,2)));
    
    % Get four box vertices
    if offCenter
        pobm = mean(pob, 1)/2;
        pbv = [-a, -b; a, -b; a, b; -a, b]/2;
        pbv = bsxfun(@plus, pbv, pobm);
        % Well center distance to the boundary
        xw =  pbv(2,1);
        yw = -pbv(2,2);
    else
        pbv = [-a, -b; a, -b; a, b; -a, b]/2;
        % Distance of well center to the boundary
        xw = a/2;
        yw = b/2;
    end
    assert(all(inpolygon(pbv(:,1), pbv(:,2), pob(:,1), pob(:,2))), ...
        ['Box vertexes outside the well reigion were detected, ', ...
        'try to reduce the box size']);
    
    % Number of angular cells
    nA  = length(bdn);
    
    % Get points on the box boundary
    pbb = zeros(nA, 2);
    sg1 = sign(pob(vxID,1:2));
    sg2 = sign(pbv);
    idx = arrayfun(@(i)find(ismember(sg1, sg2(i,:), 'rows')), (1:4)');
    pbv = pbv(idx,:);
    
    % Assign four box vertices
    pbb(vxID, :) = pbv; 
    assert(vxID(1)==1)
    for i = 1 : 4
        i1 = vxID(i);
        if i < 4
            i2 = vxID(i+1); im = (i1+1:i2-1);
        else
            i2 = 1; im = (i1+1:size(pbb,1));
        end
        p1   = pob(i1,:);
        p2   = pob(i2,:);
        pm12 = pob(im,:);
        l1m  = sqrt( sum(bsxfun(@minus, pm12, p1).^2, 2) );
        l1m  = l1m / norm(p2-p1);
        p3   = pbb(i1,:);
        p4   = pbb(i2,:);
        p34m = bsxfun(@plus, p3, l1m * (p4-p3));
        pbb(im,:) = p34m; % Assign points in box edges
    end
    
    % Get points between pob and pbb
    space = linspace(0, 1, nRT(2)+1);
    space = space(2:end-1);
    px1 = bsxfun(@plus, pob(:,1), (pbb(:,1)-pob(:,1)) * space);
    py1 = bsxfun(@plus, pob(:,2), (pbb(:,2)-pob(:,2)) * space);
    
    % Get points inside the box
    pD = getPD(rW, a, b, xw, yw, nRT, pDMult);
    [px2, py2, rij] = deal( zeros(nA, nRT(1)-1) );
    [theta, rM] = cart2pol(pbb(:,1), pbb(:,2));
    for i = 1 : size(px2,1)
        for j = 1 : size(px2,2)
            fun = @(r)computePD(r*cos(theta(i)), r*sin(theta(i)), ...
                a, b, xw, yw) - pD(j);
            rij(i,j) = bisection(fun, rW, rM(i), 1e-5);
            [px2(i,j), py2(i,j)] = pol2cart(theta(i), rij(i,j));
        end
    end
    
    % Wellbore (Casing) points
    [pxW, pyW] = pol2cart(theta, rW);
    
    % Put together, reverse for numbering
    px = [pxW, px2(:,end:-1:1), pbb(:,1), px1(:,end:-1:1), pob(:,1)];
    py = [pyW, py2(:,end:-1:1), pbb(:,2), py1(:,end:-1:1), pob(:,2)];
    
    % Return the planner points
    pxy = [px(:), py(:)];
    
    % Map back to original coordinate
    px = px(:) + origin(1);
    py = py(:) + origin(2);
    
    % Get pz by interpolation
    xx = [pob0(:,1); origin(1)];
    yy = [pob0(:,2); origin(2)];
    zz = [pob0(:,3); origin(3)];
    pz = griddata(xx, yy, zz, px, py);
    pz(1:length(bdn)) = pob0(:,3);
    
    % Return the 3D points
    p = [px, py, pz];
    p = convertTo3DPlane(p, T, R);
    
    % Store borewall and screen nodes, prepared for multi-segment well
    nA = length(bdn);
    pW = p(1:nA, :);
    wellbore.wall.radius = rW;
    wellbore.wall.coords = pW;
    wellbore.reservoirCells = (1:nA)';
    if isfield(well, 'screenRadius')
        rS  = well.screenRadius(ii);
        assert(rS < rW, ...
            'The screen radius must be less than the casing radius!')
        [pxS, pyS] = pol2cart(theta, rS);
        pxS = pxS + origin(1);
        pyS = pyS + origin(2);
        pzS = griddata(xx, yy, zz, pxS, pyS);
        pS = [pxS, pyS, pzS];
        pS = convertTo3DPlane(pS, T, R);
        wellbore.screen.radius = rS;
        wellbore.screen.coords = pS;
    end
end

% -------------------------------------------------------------------------
function pD = getPD(rW, a, b, xw, yw, nRT, pDMult)
    % Get pD inside the box
    pDM = computePD(rW, 0, a, b, xw, yw); 
    pD = linspace(pDM/pDMult, pDM, nRT(1));
    pD = pD(1:end-1);
end
