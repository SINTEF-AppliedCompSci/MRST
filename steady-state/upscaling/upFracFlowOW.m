function updata = upFracFlowOW(block, updata, varargin)
% Upscale fractional flow curves.
%
% It is assumed that all fractional flow curves are monotone.

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
    'relTolFf',    0.01 ...  % Relative to ff scale (same as absolute)
    );
opt = merge_options(opt, varargin{:});

G     = block.G;
fluid = block.fluid;

assert(isfield(fluid, 'fracFlowInv'), ...
   'The fluid structure must have the field ''fracFlowInv''. ');

pv    = block.pv;
pvTot = sum(pv);
nPointsInit = opt.nPointsInit;
nPointsMax  = opt.nPointsMax;
relTolSw    = opt.relTolSw;
relTolFf    = opt.relTolFf;

% Allocate vectors
ffWup = nan(nPointsMax, 2); % Each row is a pair of [ffup, pcOWup]

% Compute the initial set of equally spaces ff points
ffWup(1:nPointsInit, 2) = linspace(0, 1, nPointsInit)';
for i = 1:nPointsInit
   ff = ffWup(i,2)*ones(G.cells.num, 1); % Set same ff in all cells
   sw = fluid.fracFlowInv(ff); % Invert frac flow curves
   ffWup(i,1) = sum(sw.*pv) / pvTot; % Compute corresp. upscaled sat
end
sWMin = min(ffWup(1:nPointsInit, 1));
sWMax = max(ffWup(1:nPointsInit, 1));

% Continue to add points to the upscaled curve where the distance between
% two points is largest. We check distance in both directions.
nPointsExtra = nPointsMax - nPointsInit;
for i = 1:nPointsExtra+1 % Note: last iteration is just a check
   
   % Find maximum distance between two points in each direction
   [maxDiffSw, maxInxSw] = max(abs(diff(ffWup(:,1))));
   [maxDiffFf, maxInxFf] = max(abs(diff(ffWup(:,2))));
   
   % Scale the differences with the span (span of ff is 1)
   maxDiffSw = maxDiffSw/(sWMax - sWMin);
   
   if maxDiffSw < relTolSw && maxDiffFf < relTolFf
       % Both tolerences are met, and so we are done.
       if i < nPointsExtra+1
           % Remove nans at end of vectors and break out of loop
           last   = nPointsMax - (nPointsExtra-i+1);
           ffWup = ffWup(1:last,:);
       end
       break; % Break out of loop
   elseif i == nPointsExtra+1
       % Max number of iterations completed, but tolerences not met
       if mrstVerbose
           warning(['Upscaling frac flow: Tolerence not met. '...
               'sWdiff=%1.2e, ffdiff=%1.2e.'], maxDiffSw, maxDiffFf);
       end
       break; % Break out of loop
   end
   
   % Compute new pair of upscaled sw and ff. Insert where the difference
   % in either direction is largest.
   if maxDiffSw > maxDiffFf
      inx = maxInxSw;
   else
      inx = maxInxFf;
   end
   ffup = mean(ffWup([inx, inx+1],2));
   ff   = ffup*ones(G.cells.num,1); % Set same ff in all cells
   sw   = fluid.fracFlowInv(ff); % Invert fractional flow curves
   swup = sum(sw.*pv) / pvTot; % Compute corresp. upscaled sat
   
   % Insert new values at correct place to keep sorted order
   ffWup(inx+1:end,:) = [swup ffup; ffWup(inx+1:end-1,:)];
   
end

% Flip if saturations decrease
if ffWup(1,1) > ffWup(2,1)
    ffWup = flipud(ffWup);
end

% Check that upscaled values are valid
assert( (all(diff(ffWup(:,2)) > 0) || all(diff(ffWup(:,2)) < 0) ), ...
        'Upscaled fractional flow curve is not strictly monotonic');
assert( all(diff(ffWup(:,1)) > 0), ...
        'Upscaled water saturation is not strictly increasing');

% Add upscaled data to structure
updata.ffW = ffWup;

end
