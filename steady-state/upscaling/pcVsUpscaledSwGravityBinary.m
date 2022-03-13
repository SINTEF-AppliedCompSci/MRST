function [pcFun, swUMin, swUMax, pcowVec] = pcVsUpscaledSwGravityBinary(...
   G, rock, fluid, varargin)
% Create fine scale capillary pressure pc as a function of upscaled sW.
% 
% REQUIRED PARAMETERS:
%   G        - Grid structure.
%   rock     - Rock structure.
%   fluid    - Fluid structure. Must have a field called 'pcOWInv'.
% 
% RETURNS:
%   pcFun    - function handle to a function taking upscaled sW as input, 
%              and returning the corresponding capillary pressure.
%   swUMin   - minimum upscaled water saturation
%   swUMax   - maximum upscaled water saturation

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

% opt = struct(...
%    'nPointsMin',   30, ...
%    'nPointsMax',  550, ... % Overrules difftols (breaks after nPointsMax)
%    'diffTolSw',  0.001, ...
%    'diffTolPc',  0.001 ... % Relative to pc scale
%    );
% opt = merge_options(opt, varargin{:});

assert(isfield(fluid, 'pcOWInv'), ...
   'The fluid structure must have the field ''pcOWInv''. ');

% ----------------------------------------------------------

% Functions
pcOW = fluid.pcOW;
pcOWInv = fluid.pcOWInv;

% Pore volume
pv = G.cells.volumes.*rock.poro;

% Gravity force value for each cell in the grid
g = 9.8066; % HARDCODED
if isfield(fluid, 'rhoO')
   rhoO = fluid.rhoO;
elseif isfield(fluid, 'rhoOS')
   rhoO = fluid.rhoOS;
else
   error('No oil denisty field in fluid');
end
if isfield(fluid, 'rhoW')
   rhoW = fluid.rhoW;
elseif isfield(fluid, 'rhoWS')
   rhoW = fluid.rhoWS;
else
   error('No water denisty field in fluid');
end
dRho = rhoO - rhoW;

% Compute an estimate of the centroid of the grid block. This will be the
% correct centroid is the grid is Cartesian, but otherwise, it may be off.
zCent = mean([max(G.cells.centroids(:,3)), min(G.cells.centroids(:,3))]);
%zCent = max(G.cells.centroids(:,3));
%zCent = min(G.cells.centroids(:,3));

% Set height relative to the zCent. Thus the returned xvec is the capillary
% pressure at the height zCent.
zi   = zCent - G.cells.centroids(:,3);

grav = dRho.*g.*zi;

% Assume we are given swir and sor, and also assume they are the same for
% all rock types
assert(isfield(fluid, 'swir') && isfield(fluid, 'sor'), ...
   'Need fluid fields ''swir'' and ''sor''.');
assert(all(fluid.swir==fluid.swir(1)) && all(fluid.sor==fluid.sor(1)), ...
   'Different saturation limits not supported.');
swlim = [fluid.swir, 1-fluid.sor];

% Compute relationship between upscaled saturation and the capillary
% pressure (at a some datum)
[avgSW, pcowVec] = preComputeSw(pcOW, pcOWInv, pv, ...
   grav, swlim);

% For a given upscaled water saturation, return capillary pressure
pcFun = @(sWU) interp1(avgSW, pcowVec, sWU);

pcowVec = [avgSW pcowVec];

swUMin = min(avgSW);
swUMax = max(avgSW);

end



function [avgSW, xVec] = preComputeSw(pcOW, pcOWInv, pv, ...
   grav, swlim, varargin)
% We precompute the relation between the mean saturation and the cell
% saturation. Then, we can later do lookup during the time integration.

opt = struct(...
   'method',    'binary',  ...
   'tol',       1e-10, ...
   'npoints',   200 ...
   );
opt = merge_options(opt, varargin{:});

% Number of points
P = opt.npoints; % Total number of points (including initial points)

% Allocation
avgSW = nan(P, 1);
xVec  = nan(P,1);
sW    = nan(numel(pv), 1);

% Compute left outer point, sw=swir
sWmean = swlim(1);
sW(:)  = sWmean;
xVal   = 0;
[~, xVal] = findSwNewton(sW, xVal, sWmean, ...
   pcOW, pcOWInv, pv, grav, swlim, 'tol', opt.tol);
xVec(1)  = xVal;
avgSW(1) = sWmean;

% Compute right outer point, sw=1-sor
sWmean = swlim(2);
sW(:)  = sWmean;
xVal   = 0;
[~, xVal] = findSwNewton(sW, xVal, sWmean, ...
   pcOW, pcOWInv, pv, grav, swlim, 'tol', opt.tol);
xVec(2)  = xVal;
avgSW(2) = sWmean;

% Compute the span
sSpan = abs(diff(avgSW(1:2)));
xSpan = abs(diff(xVec(1:2)));

% Compute remaining points in between the outer points
for p = 3:P
   
   % We insert a new point at the largest inteval among the points.
   % The relative distance in both sw and x is checked.
   [sDiff, sInx] = max(abs(diff(avgSW)));
   [xDiff, xInx] = max(abs(diff(xVec)));
   if sDiff/sSpan >= xDiff/xSpan
      inx = sInx;
   else
      inx = xInx;
   end
   sWmean  = mean(avgSW(inx:inx+1));
   
   % Find saturation distribution
   xilim = xVec(inx:inx+1);
   [~, xVal] = findSwBinary(sWmean, pcOWInv, pv, grav, ...
      xilim, 'tol', opt.tol);
   
   % Insert computed value into vectors and shift remaining part of vector
   xVec(inx+1:p)     = [xVal;      xVec(inx+1:p-1)];
   avgSW(inx+1:p)    = [sWmean;   avgSW(inx+1:p-1)];
end

end


function [sW, xi, meta] = findSwBinary(sWmean, pcOWInv, pv, grav, ...
   xilim, varargin)
% The average saturation is given, and then the saturation in each cell is
% computed, as well as the capillary pressure.
opt = struct(...
   'tol',  1e-10  ...
   );
opt = merge_options(opt, varargin{:});

% Test that the limits are given in the right order
sWlim(1) = sum(pv.*pcOWInv(xilim(1) - grav))/sum(pv);
sWlim(2) = sum(pv.*pcOWInv(xilim(2) - grav))/sum(pv);
if sWlim(1) > sWlim(2)
   sWlim = sWlim([2,1]);
   xilim = xilim([2,1]);
end
diffs = sWlim - sWmean;

xi  = 0;
itr = 0;
tol = opt.tol;
maxItr = 50;
conv = false;
while itr < maxItr
   
   itr = itr + 1;
   sW  = pcOWInv(xi - grav);
   diff = sum(pv.*sW)/sum(pv) - sWmean;
   
   % Check for convergence
   err = abs(diff);
   if err < tol
      conv = true;
      break
   end
   
   if diff > 0
      % Saturation is too large
      xilim(2) = xi; % New upper boundary
      diffs(2) = diff;
   else
      % Saturation is too small
      xilim(1) = xi; % New lower boundary
      diffs(1) = diff;
   end
   
   % Make a new guess half way in between
   xi = xilim(1) + 0.5*(xilim(2)-xilim(1));
   
end

% Check convergence
if ~conv
   warning(['The solver did not converge for '...
      'sWmean=%1.4f.\nThe convergence error was %1.2e.'], ...
      sWmean, err );
else
%    fprintf(['Converged for sWmean=%1.2f after %d itr.\n' ...
%       'The convergence error was %1.2e.\n'], sWmean, itr, err);
end

meta.err = err;
meta.itr = itr;

end


function [sW, xVal, meta] = findSwNewton(sW, xVal, sWmean, ...
   pcOW, pcOWInv, pv, grav, swlim, varargin)
% The average saturation is given, and then the saturation in each cell is
% computed, as well as the capillary pressure.
opt = struct(...
   'tol',  1e-10  ...
   );
opt = merge_options(opt, varargin{:});

% Settings
tol    = opt.tol;
maxItr = 20;
dsMax = 0.05;

% We have n unknowns in sW, and 1 unknown in pcval. To find these, we solve
% n+1 nonlinear equations using a very simple Newton loop. This is done
% using ADI.

conv = false;
itr  = 0;
while ~conv && itr < maxItr
   
   % Iteration counter
   itr    = itr + 1;
   
   % Create ADI variables
   [sW, xVal] = initVariablesADI(sW, xVal);
   
   % Equation for pcow and gravity
   % The domain is divided into two, where the equation is used on the one
   % section, and the inverse problem on the other.
   soff = 0.05;
   inx  = sW<swlim(1)+soff | sW>swlim(2)-soff;
   invVal = pcOWInv(xVal - grav);
   eqs{1} = (pcOW(sW) + grav - xVal).*(1/barsa); % scaling
   eqs{1}(inx) = sW(inx) - invVal(inx); % edges
   
   % Mean saturation equation
   eqs{2} = sWmean - sum(pv.*sW) / sum(pv);
   
   % Check for convergence
   resVals = [max(abs(eqs{1}.val)), max(abs(eqs{2}.val))];
   conv    = max(resVals) < tol;
   if conv
      sW     = sW.val;
      xVal  = xVal.val;
      break
   end
   
   % Solve
   dx     = SolveEqsADI(eqs, []);
   
   if ~iscell(dx)
      fprintf('Fail at itr=%d for sWmean=%1.2f\n', itr, sWmean);
   end
   
   % Get increments
   ds     = dx{1};
   dp     = dx{2};
   
   % Limit increments
   ds = sign(ds).*min(abs(ds), dsMax);

   % Apply increments
   sW     = sW.val + ds;
   xVal  = xVal.val + dp;
   
   % Enforce limits
   sW     = min(swlim(2), max(swlim(1), sW));
   
end

% Check convergence
if ~conv
   warning(['The nonlinear solver did not converge for '...
      'sWmean=%1.2f.\nThe convergence error was %1.2e.'], ...
      sWmean, max(resVals) );
else
%    fprintf(['Converged for sWmean=%1.2f after %d itr.\n' ...
%       'The convergence error was %1.2e. '], sWmean, itr, max(conVals));
%    fprintf('Also, sum(inx)=%d\n', sum(inx));
end

meta.converged  = conv;
meta.iterations = itr;
meta.residuals  = resVals;

end
