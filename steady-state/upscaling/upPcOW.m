function [updata, report] = upPcOW(block, updata, varargin)
% Upscale capillary pressure curves.
%
% It is assumed that all capillary pressure curves are monotone.
% 
% The upscaling is done my looping over different values of the capillery
% pressure. For each value, this value is set in all cells, and then the
% capillary pressure function is inverted to find the saturation
% distribution. The saturation is then upscaled, and we have a pair of
% saturation/pcOW.

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

opt = struct(...
    'nPointsInit', 20, ...   % Number of initial points
    'nPointsMax',  50, ...   % Maximum number of points
    'relTolSw',    0.01, ... % Relative to saturation scale
    'relTolPc',    0.01, ...  % Relative to pc scale
    'gravity',     'none', ...
    'fixjumps',    true ...
    );
opt = merge_options(opt, varargin{:});

wantReport = nargout > 1;
timeStart = tic;

G     = block.G;
fluid = block.fluid;

if ~isfield(fluid, 'pcOW')
	updata.pcOW = [];
    return
end

assert(isfield(fluid, 'pcOWInv'), ...
   'The fluid structure must have the field ''pcOWInv''. ');

pv    = block.pv;
pvTot = sum(pv);
nPointsInit = opt.nPointsInit;
nPointsMax  = opt.nPointsMax;
relTolSw    = opt.relTolSw;
relTolPc    = opt.relTolPc;

% Handle gravity
useGravity  = false;
grav        = 0;
if any(strcmpi(opt.gravity, {'centroid', 'bottom'}))
    useGravity = true;
    
    if isfield(fluid, 'rhoO')
        rhoO = fluid.rhoO;
    else
        rhoO = fluid.rhoOS;
    end
    if isfield(fluid, 'rhoW')
        rhoW = fluid.rhoW;
    else
        rhoW = fluid.rhoWS;
    end
    dRho = rhoW - rhoO;
    g    = 9.8066; % HARDCODED
    
    if strcmpi(opt.gravity, 'centroid')
        % Compute an estimate of the centroid of the grid block. This will
        % be the correct centroid is the grid is Cartesian, but otherwise,
        % it may be off.
        zCent = mean([max(G.cells.centroids(:,3)), ...
                      min(G.cells.centroids(:,3))]);

        % Set height relative to the zCent. Thus the returned xvec is the
        % capillary pressure at the height zCent.
        zi   = zCent - G.cells.centroids(:,3);
    elseif strcmpi(opt.gravity, 'bottom')
        zi   = max(G.cells.centroids(:,3)) - G.cells.centroids(:,3);
    end
    
    grav = -dRho.*g.*zi; % NOTE: Not sure about sign here
end

% Find minimum and maximum pc in block
pc0   = fluid.pcOW( zeros(G.cells.num,1) ) + grav;
pc1   = fluid.pcOW( ones(G.cells.num,1) ) + grav;
pcMin = min(min(pc0), min(pc1));
pcMax = max(max(pc0), max(pc1));

% Allocate vectors
pcOWup = nan(nPointsMax, 2); % Each row is a pair of [sWup, pcOWup]

% Compute the initial set of equally spaces pc points
pcOWup(1:nPointsInit, 2) = linspace(pcMax, pcMin, nPointsInit)';
for i = 1:nPointsInit
   pc = pcOWup(i,2)*ones(G.cells.num, 1); % Set same pc in all cells
   sw = fluid.pcOWInv(pc - grav); % Invert capillary pressure curves
   pcOWup(i,1) = sum(sw.*pv) / pvTot; % Compute corresp. upscaled sat
end
sWMin = min(pcOWup(1:nPointsInit, 1));
sWMax = max(pcOWup(1:nPointsInit, 1));

% Continue to add points to the upscaled curve where the distance between
% two points is largest. We check distance in both directions.
nPointsExtra = nPointsMax - nPointsInit;
for i = 1:nPointsExtra+1 % Note: last iteration is just a check
   
   % Find maximum distance between two points in each direction
   [maxDiffSw, maxInxSw] = max(abs(diff(pcOWup(:,1))));
   [maxDiffPc, maxInxPc] = max(abs(diff(pcOWup(:,2))));
   
   % Scale the differences with the span
   maxDiffSw = maxDiffSw/(sWMax - sWMin);
   maxDiffPc = maxDiffPc/(pcMax - pcMin);
   
   if maxDiffSw < relTolSw && maxDiffPc < relTolPc
       % Both tolerences are met, and so we are done.
       if i < nPointsExtra+1
           % Remove nans at end of vectors and break out of loop
           last   = nPointsMax - (nPointsExtra-i+1);
           pcOWup = pcOWup(1:last,:);
       end
       break; % Break out of loop
   elseif i == nPointsExtra+1
       % Max number of iterations completed, but tolerences not met
       if mrstVerbose
           warning(['Upscaling pcOW: Tolerence not met. '...
               'sWdiff=%1.2e, pcdiff=%1.2e.'], maxDiffSw, maxDiffPc);
       end
       break; % Break out of loop
   end
   
   % Compute new pair of upscaled sw and pcow. Insert where the difference
   % in either direction is largest.
   if maxDiffSw > maxDiffPc
      inx = maxInxSw;
   else
      inx = maxInxPc;
   end
   pcup = mean(pcOWup([inx, inx+1],2));
   pc   = pcup*ones(G.cells.num,1); % Set same pc in all cells
   sw   = fluid.pcOWInv(pc - grav); % Invert capillary pressure curves
   swup = sum(sw.*pv) / pvTot; % Compute corresp. upscaled sat
   
   % Insert new values at correct place to keep sorted order
   pcOWup(inx+1:end,:) = [swup pcup; pcOWup(inx+1:end-1,:)];
   
end

% Flip if saturations decrease
% We use a bit silly method to determine if it is deacreasing or not. We
% used to check only the first two values, but they can in some cases be
% zero, and so it is not enough to check.
if sum(diff(pcOWup(:,1))<0) > size(pcOWup,1)/2
    pcOWup = flipud(pcOWup);
end

if opt.fixjumps
    % Remove consequitive points with same upscaled saturation
    jumpinx  = [true; diff(pcOWup(:,1))~=0];
    pcOWup   = pcOWup(jumpinx, :);
    njumpfix = sum(~jumpinx);
    if njumpfix>0
        %warning('Removed %d points in pcOW due to jumps.', njumpfix);
    end
end

% Check that upscaled values are valid
assert( (all(diff(pcOWup(:,2)) > 0) || all(diff(pcOWup(:,2)) < 0) ), ...
        'Upscaled capillary pressure curve is not strictly monotonic');
assert( all(diff(pcOWup(:,1)) > 0), ...
        'Upscaled water saturation is not strictly increasing');

% Add upscaled data to structure
updata.pcOW = pcOWup;

totalTime = toc(timeStart);
if wantReport
    report.swreldiff = maxDiffSw;
    report.pcreldiff = maxDiffPc;
    report.gravity   = useGravity;
    report.time      = totalTime;
    if opt.fixjumps
        report.numjumpfix = njumpfix;
    end
end

end
