classdef WellPositionControl < handle
    % Class for position control for single well
    properties
       w                                % scalar well structure
       nPoints              = 2;        % number of interpolation-points
       controlPoints                    % 3D interpolation points
       activePoints                     % which points are alowed to move
       linkedPoints                     % option for linking points to other instances
       parameters                       % set up by contructor
       perturbationSize                 % x,y and z axis of perturbaion ellipsoid
       maxUpdateWell        = [inf, inf, inf];
       maxUpdatePoint       = [inf, inf, inf];
       maxRadianChange      = (1/3)*2*pi/360;     
       minimalLength        = 5*meter;
       regionIx
       regions
       G
       monotoneZ            = true;     % option to enforce monotonicity in z-coord
       boxes                = {};       % bounding box of each regian (used for scaling)
       regionTollerance     =  1*meter; % projection tollerance (+/- inside/outside)
       interpolationMethod  = 'pchip';
       trajectoryTollerance = .5*meter; % traj tollerance (max diff between curve and piecewise linear approximation)
       its                  = struct('nOuter', 1, ...
                                     'nInner', 1);
    end
    methods 
        function p = WellPositionControl(G, varargin)
            p = merge_options(p, varargin{:});
            if ~isempty(p.w)
                % setup from exising well
                p = setupFromExistingWell(p, G);
            else
                if ~isempty(p.controlPoints)
                    p.nPoints = size(p.controlPoints,1);
                end
                if isempty(p.activePoints)
                    p.activePoints = true(p.nPoints, 1);
                end
                if isempty(p.regionIx)
                    p.regionIx = ones(p.nPoints, 1);
                end
                if isempty(p.regions)
                    p.regions{1} = boundaryFaces(G);
                end
                if numel(p.regionTollerance) == 1
                    p.regionTollerance = repmat(p.regionTollerance, [p.nPoints, 1]);
                end
            end
            p = setupParameters(p);
            if isempty(p.boxes)
                p.boxes = cell(1, p.nPoints);
                for k = 1:p.nPoints
                    p.boxes{k} = getRegionBoundingBox(p,G, k);
                end
            end
            p.G = G;
        end
        
        function control = param2control(p, p1)
            % realative controlPoints in box
            if nargin < 2
                p1 = p.controlPoints;
            end
            control = cell(1, p.nPoints);
            assert(p.nPoints*3 == p.parameters.nParam, 'for now');
            for k = 1:p.nPoints
                %[mn, mx] = deal(p.boxes{k}(1,:), p.boxes{k}(2,:));
                %ck = (p1(k,:)-mn)./(mx-mn);
                ck = p1(k,:);
                control{k} = ck(:);
            end
            control = vertcat(control{:});
        end
        
        function p = control2param(p, control)
            %if nargin < 3
            %    force = false;
            %end
            assert(numel(control) == p.parameters.nParam);
            ix = 0;
            newPoints = nan(size(p.controlPoints));
            for k = 1:p.nPoints
                %[mn, mx] = deal(p.boxes{k}(1,:), p.boxes{k}(2,:));
                newPoints(k,:) = control(ix + (1:3))';%.*(mx-mn)+mn;
                ix = ix +3;
            end
            if true
                p.controlPoints = newPoints;
                p = setupParameters(p);
            else
                
                p = setupParameters(p);
            end
        end
        
        function t = getTrajectory(p)
            % sample points along trajectory that confine with tollerance
            f = p.getInterpolationFun;
            t = samplePiecewiseLinear(f, (0:1)', p.trajectoryTollerance);
        end
        
        function f = getInterpolationFun(p)
            dp     = diff(p.controlPoints);
            cumlen = [0; cumsum( sqrt( dot(dp, dp, 2) ) )];
            f = @(s)interp1(cumlen/cumlen(end), p.controlPoints, s, ...
                            p.interpolationMethod, 'extrap');
        end
        
        function du = getProjectedUpdate(p, u, du, isScaled)
            dp = nan(size(p.controlPoints));
            ix = 0;
            for k = 1:p.nPoints
                if isScaled
                    [mn, mx] = deal(p.boxes{k}(1,:), p.boxes{k}(2,:));
                else
                    [mn, mx] = deal(0,1);
                end
                dp(k,:) = du(ix + (1:3))'.*(mx-mn);
                ix = ix +3;
            end
            dp = getProjectedPointUpdate(p, dp);
            up = param2control(p, p.controlPoints + dp);
            ix = 0;
            for k = 1:p.nPoints
                if isScaled
                    [mn, mx] = deal(p.boxes{k}(1,:), p.boxes{k}(2,:));
                else
                    [mn, mx] = deal(0,1);
                end
                du(ix + (1:3)') = dp(k,:)./(mx-mn);
                ix = ix +3;
            end
            %du = up -u;
        end
            
        function dpp = getProjectedPointUpdate(p, dp)
            if norm(dp) < sqrt(eps)*max(norm(p.controlPoints), 1)
                %warning('Requested zero update');
                return
            end
            % Outer:
            % 1 Move all active points in direction min(sum(dp(active), maxUpdateWell))
            %    a) if some points outside, move all according worst
            %       offender
            %    b) enforce monotonicity  (one sweep from heel to toe)
            %    c) enforce max curvature (one sweep from heel to toe)
            %    record update v_outer
            %
            % Inner
            %    For each point (from toe to heel) consider
            %       vk = dp(k,:)-v_outer(k, :)
            %    Move point in direcction min(vk, maxUpdatePoint))
            %    a) if point is outside, move to boundary
            %    b) enforce monotonicity  (if not first)
            %    c) enforce max curvature (if not first/last)
            %
            %   Iterate if nIts > 1
            % 1 max update
            
            % MOVE WHOLE SET OF ACTIVE POINTS
            pnt = p.controlPoints;
            for outer = 1:p.its.nOuter
                aix  = p.activePoints & ~p.linkedPoints;
                vec  =  absMin(mean(dp(aix,:)), p.maxUpdateWell);
                if all(aix) % we can move freely
                    upd = aix*vec;
                else
                    % do curvature update:
                    tmp = projectRadian(pnt + vec, p.maxRadianChange, pnt, ~aix);
                    upd = pnt - tmp;
                end
                % check regions
                for k = 1:p.nPoints
                    if aix(k)
                        reg = p.regions{p.regionIx(k)};
                        d = getDistanceToBoundary(p.G, pnt(k,:), upd(k,:), reg);
                        if d - p.regionTollerance <= norm(upd(k,:))
                            tmp = getDirectionToClosestFace(p.G, pnt(k,:)+upd(k,:), reg);
                            upd(k,:) = upd(k,:) + tmp*(1 + p.regionTollerance(k)/norm(tmp));
                        end
                    end
                end
                pnt = pnt + upd;
                pnt = enforceMonotoneZ(pnt, aix, p.monotoneZ);
                pnt = enforceMinLength(pnt, aix, p.minimalLength);
                for inner = 1:p.its.nInner
                    % MOVE INDIVIDUAL POINTS
                    upd = dp - (pnt - p.controlPoints);
                    for k = 1:p.nPoints
                        upd(k,:) = absMin(upd(k,:), p.maxUpdatePoint);
                    end
                    for k = 1:p.nPoints
                        if aix(k)     
                            % curvature
                            if k > 1 && k < p.nPoints
                                ix = k-1:k+1;
                                aix_cur = aix(ix) & logical([0 1 0]');
                                tmp = projectRadian(pnt(ix,:) + upd(ix,:), ...
                                    p.maxRadianChange, pnt(ix,:), ~aix_cur);
                                upd(k,:) = tmp(2,:)-pnt(k,:);
                            end
                            % regions
                            reg = p.regions{p.regionIx(k)};
                            d = getDistanceToBoundary(p.G, pnt(k,:), upd(k,:), reg);
                            if d - p.regionTollerance <= norm(upd(k,:))
                                tmp = getDirectionToClosestFace(p.G, pnt(k,:)+upd(k,:), reg);
                                upd(k,:) = upd(k,:) + tmp*(1 + p.regionTollerance(k)/norm(tmp));
                            end
                        end
                    end
                    
                    pnt = pnt + upd;
                    pnt = enforceMonotoneZ(pnt, aix, p.monotoneZ);
                    pnt = enforceMinLength(pnt, aix, p.minimalLength);
                end
            end
            dpp = pnt - p.controlPoints;
        end
        
      
    end
end

% -------------------------------------------------------------------------

function p = setupFromExistingWell(p, G)
% assert(isempty(p.nPoints) || isempty(p.controlPoints))
if isempty(p.controlPoints)
    assert(p.nPoints >= 2);
    % interpolate along trajectory
    if ~isfield(p.w, 'traj') % we need to approximate traj from connection cells
        tmp = addTrajectories(p.w, G, 10);
        % check badness ??
    end
    % position controlPoints equally along traj
    traj = tmp.trajectory;
    seg  = diff(traj);
    cumlen = [0; cumsum( sqrt( dot(seg,seg,2) ) )];
    p.controlPoints = interp1(cumlen/cumlen(end), traj, (0:p.nPoints-1)/(p.nPoints-1), 'linear', 'extrap');
    % update well-trajectory according to controlPoints and interp-method
    p.w.traj = p.getTrajectory;
end
% assert(p.nPoints == size(p.controlPoints,1))
if isempty(p.activePoints)
    p.activePoints = true(p.nPoints, 1);
end
if isempty(p.linkedPoints)
    p.linkedPoints = false(p.nPoints, 1);
end
if ~isempty(p.regionIx)
    assert(max(p.regionIx) <= numel(p.regions));
else
    p.regionIx = ones(p.nPoints, 1); % assume within first region
end
if isempty(p.regions) % assume whole grid
    p.regions = {boundaryFaces(G)};
end
if numel(p.regionTollerance) == 1
    p.regionTollerance = repmat(p.regionTollerance, [p.nPoints, 1]);
end
p.maxUpdatePoint = min(p.maxUpdatePoint, p.maxUpdateWell);
% remove cells from well (it's a template and its cells should not be used)
p.w.cells = [];
end

% -------------------------------------------------------------------------

function p = setupParameters(p)
assert(~isempty(p.perturbationSize))
r = p.perturbationSize;
dp     = diff(p.controlPoints);
cumlen = [0; cumsum( sqrt( dot(dp, dp, 2) ) )];
[pointIx, vec, lengthDerivative] = deal(cell(p.nPoints, 1));
for k = 1:p.nPoints
    if p.activePoints(k)
        % get trajecroty tangent at point
        f  = p.getInterpolationFun();
        sk = cumlen(k)/cumlen(end);
        Pv = diff(f(sk + sqrt(eps)*[-1, 1]'));
        Pv = Pv/norm(Pv);
        
        % normal plane: Pv*x' = 0
        % ellipsoid: sum( (x./r).^2 ) = 1
        % change of variables:
        % u -> x/rx, v -> y/ry, w -> z/rz
        
        % plane in (u,v,w)-system: m*x' = 0
        m = Pv.*r;
        m = m/norm(m);
        if norm(m(1:2)) < sqrt(eps) % horizontal plane
            [e1, e2] = deal([1 0 0], [0 1 0]);
        else
            e1 = [m(2), -m(1), 0]/norm(m(1:2));
            e2 = cross(m, e1);
        end
        % change back to (x,y,z)
        [v1, v2] = deal(e1.*r, e2.*r);
        v = [v1;v2];
        isEnd = k==1 || k==p.nPoints;
        lengthDerivative{k} = zeros(2+isEnd,1);
        if isEnd
            v3 = m.*r;
            v  = [v;v3]; %#ok
            ld = sign( dot(Pv, v3) );
            if k==1
                ld = -ld;
            end
            lengthDerivative{k}(3) = ld;
        else
            
        end
        pointIx{k} = k*ones(2+isEnd,1);
        vec{k}     = v;
    end
end
pointIx = vertcat(pointIx{:});
p.parameters = struct('pointIx', pointIx, 'perturbation', vertcat(vec{:}), ...
                      'lengthDerivative', vertcat(lengthDerivative{:}), ...
                      'nParam', 3*max(pointIx));
end

% -------------------------------------------------------------------------

function box = getRegionBoundingBox(p, G, k)
fc = G.faces.centroids(p.regions{p.regionIx(k)}, :);
box = [min(fc)+p.regionTollerance(k); max(fc)-p.regionTollerance(k)];
end
% -------------------------------------------------------------------------
function a = absMin(a, b)
    if ~isempty(b)
        a = sign(a).*min(abs(a), abs(b));
    end
end

% -------------------------------------------------------------------------
function pnt = enforceMonotoneZ(pnt, ix, do)
if do
    for k = 2:size(pnt,1)
        if ix(k)
            pnt(k,3) = max(pnt(k,3), pnt(k-1,3));
        end
    end
end
end
% -------------------------------------------------------------------------
function pnt = enforceMinLength(pnt, ix, minLength)
    np = size(pnt,1);
    dp = diff(pnt);
    l  = sqrt(dot(dp,dp,2));
    if sum(l) < minLength
        assert(ix(end))
        vec = sum(dp,1);
        if norm(vec) < sqrt(eps)
            vec = [0 0 1];
        else
            vec = vec/norm(vec);
        end
        pnt = pnt(1,:) + minLength*((0:np-1)'*vec);
    elseif any(l) < sqrt(eps)
        for k = 2:np-1
            if l(k-1) < sqrt(eps) && ix(k)
                pnt(k,:) = .5*(pnt(k-1,:)+pnt(min(k+1, np)));
            end
        end
    end
end


%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

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

