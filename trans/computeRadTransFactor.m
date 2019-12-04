function ft = computeRadTransFactor(G, pW, skin, varargin)
% Compute the radial half transmissibility factor (ft) for the 2D radial 
% grid (halfTrans = perm .* ft). The computation assumes the steady-state 
% flow, and the 'transmissibility center' is obtained by integral of the 
% pressure within the area of cell
%
% SYNOPSIS:
%    ft = computeRadTransFactor(G, pW, skin)
%    ft = computeRadTransFactor(G, pW, skin, 'nodeCoords', nodeCoords)
%
% PARAMETERS:
%    G     - Radial grid structure, typically built by `buildRadialGrid`. 
%            The numbering of cells and nodes obey the the logical
%            numbering (angular cycle fastest, then radial).
%            G should contain the field: 'radDims': [nA, nR] or 
%            [nA, nR(1), nR(2)]. For the second case, the total raidal 
%            dimension is nR(1)+nR(2), but only cells with r-indices within 
%            1 - nR(1) are involved in the calculations.
%    pW    - 2D coordinate of the well point
%    skin  - Skin factor of the well
%
% KEYWORD ARGUMENTS:
%   'nodeCoords'  - Provided 2D coordinates of grid nodes. The G will be
%                   updated after assigning 'G.nodes.coords' = nodeCoords.
%
% RETURNS:
%    ft   - Radial half transmissibility factor, corresponding to
%           'G.cells.faces'. ft of cells with r-indices within 
%            nR(1)+1 - nR(1)+nR(2) are 'nan'.
%
% EXAMPLE:
%   [nA, nR, rW, rM] = deal(40, 10, 2, 10);
%   th = linspace(0, 2*pi, nA+1); th = th(1:end-1);
%   r = logspace(log10(rW), log10(rM), nR+1);
%   [R, TH] = meshgrid(r, th);
%   [px, py] = pol2cart(TH(:), R(:));
%   p = [px(:), py(:)];
%   G = buildRadialGrid(p, nA, nR);
%   ft = computeRadTransFactor(G, [0, 0], 0);
%
% SEE ALSO:
%   `computeTrans`

    opt = struct('nodeCoords', []);
    opt = merge_options(opt, varargin{:});
    
    if ~isempty(opt.nodeCoords)
        p = opt.nodeCoords;
    else
        p = G.nodes.coords;
    end
    G.nodes.coords = p;
    G = computeGeometry(G);
    [nA, nRT]  = deal(G.radDims(1), G.radDims(2:end));
    ft = nan(nA*sum(nRT)*4, 1); 
    pA = pW;
    % nodes:
    % B1 : R- & A+, 1 & 2
    % C1 : R- & A-, 1 & 4
    % C2 : R+ & A-, 3 & 4
    % B2 : R+ & A+, 2 & 3
    N = [1, 2; 1, 4; 3, 4; 2, 3];
    for c = 1 : nA*nRT(1)
        fPos = G.cells.facePos(c) : G.cells.facePos(c+1)-1;
        f = G.cells.faces(fPos);
        n = arrayfun(@(f)G.faces.nodes(G.faces.nodePos(f): ...
            G.faces.nodePos(f+1)-1), f, 'UniformOutput', false);
        nBC = zeros(4,1);
        for i = 1 : 4
            nBC(i) = intersect(n{N(i,1)}, n{N(i,2)});
        end
        pB1 = p(nBC(1), :);
        pC1 = p(nBC(2), :);
        pC2 = p(nBC(3), :);
        pB2 = p(nBC(4), :);
        [IR1, IA1, S1,   ~, b1, c1] = computeI(pA, pB1, pC1);
        [IR2, IA2, S2, thA, b2, c2] = computeI(pA, pB2, pC2);
        S0   = S2 - S1;
        r0   = exp( S2/S0*IR2 - S1/S0*IR1 );
        dth0 = S2/S0*IA2 - S1/S0*IA1;
        r_ = norm(G.faces.centroids(f(1), :)-pA, 2);
        r  = norm(G.faces.centroids(f(3), :)-pA, 2);
        if c <= nA
            % Equivalent radius
            r_ = r_ * exp(-skin);
        end
        ft_r_  = thA / log(r0/r_);          % r-
        ft_r   = thA / log(r/r0);           % r+
        ft_a_B = log(b2/b1) / (thA - dth0); % th B
        ft_a_C = log(c2/c1) / dth0;         % th C
        
        % face type:
        % 1 - Radial -
        % 2 - Angular B1B2
        % 3 - Radial +
        % 4 - Angular C1C2
        ft(fPos)  = [ft_r_;    ft_a_B;    ft_r;      ft_a_C];
        if any(ft(fPos) < 0)
            1;
        end
    end
end

function [IR, IA, S, thA, b, c] = computeI(pA, pB, pC)
    fAngle = @(v1,v2)acos(dot(v1,v2)/norm(v1,2)/norm(v2,2));
    
    thA = fAngle(pB-pA, pC-pA);
    thB = fAngle(pA-pB, pC-pB);
    thC = fAngle(pB-pC, pA-pC);

    c = norm(pA-pB,2);
    b = norm(pA-pC,2);
    a = norm(pB-pC,2);

    IR = c/a*sin(thB)*thA...
        + b/a*cos(thC)*log(b)...
        + c/a*cos(thB)*log(c)...
        -1.5;

    IA = b/a*sin(thC)*...
        ( log(sin(thC)/sin(thB)) + (pi/2 - thC)*cot(thC) ...
        - (pi/2 - thB)*cot(thB) ) + pi/2 - thB;

    S = tri_area(pA, pB, pC);
end
