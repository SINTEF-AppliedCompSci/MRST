function [ffFun, swUMin, swUMax] = fracFlowVsUpscaledSw(G, rock, fluid, ...
   varargin)
% Create fine scale fractional flow as a function of upscaled sW.
% 
% REQUIRED PARAMETERS:
%   G        - Grid structure.
%   rock     - Rock structure.
%   fluid    - Fluid structure. Must have a field called 'fracFlowInv'.
% 
% RETURNS:
%   ffFun    - function handle to a function taking upscaled sW as input, 
%              and returning the corresponding fractional flow.
%   swUMin   - minimum upscaled water saturation
%   swUMax   - maximum upscaled water saturation

opt = struct(...
   'nPointsMin',   30, ...
   'nPointsMax',   50, ... % Overrules difftols
   'diffTolSw',  0.05, ...
   'diffTolFf',  0.01 ...
   );
opt = merge_options(opt, varargin{:});

assert(isfield(fluid, 'fracFlowInv'), ...
   'The fluid structure must have the field ''fracFlowInv''.');

% Options
nPointsMin = opt.nPointsMin; % initial ff data points
nPointsMax = opt.nPointsMax; % max number of extra data points
diffTolSw  = opt.diffTolSw; % diff tolerence for upscaled sw
diffTolFf  = opt.diffTolFf; % diff tolerence for factional flow

% Allocate vectors
ffVec = nan(nPointsMax, 1);
swU   = nan(nPointsMax, 1); % Upscaled water sat

% Compute the initial set of equally spaces pc points
ffVec(1:nPointsMin) = linspace(0, 1, nPointsMin)';
pv    = rock.poro.*G.cells.volumes;
pvTot = sum(pv);
for i = 1:nPointsMin
   ff = ffVec(i)*ones(G.cells.num, 1);
   sw = fluid.fracFlowInv(ff); % inverse
   swU(i) = sum(sw.*pv) / pvTot;
end
swUMin = min(swU);
swUMax = max(swU);

% Continue to add points to the curve where the distance between two points
% is large.
nPointsExtra = nPointsMax - nPointsMin;
for i = 1:nPointsExtra+1 % Note: last iteration is just a check
   
   % Check distance between points
   [maxDiffSw, maxInxSw] = max(abs(diff(swU)));
   [maxDiffPc, maxInxFf] = max(abs(diff(ffVec)));
   if maxDiffSw < diffTolSw && maxDiffPc < diffTolFf
      if i == nPointsExtra+1
         break; % ok. Reached tolerence in last iteration
      else
         % Remove nans at end of vectors and break out of loop.
         last  = length(ffVec) - (nPointsExtra-i+1);
         ffVec = ffVec(1:last);
         swU  = swU(1:last);
      end
      break;
   elseif i == nPointsExtra+1
      warning('Max number of points reached, but tolerence not met.');
      break;
   end
   
   % Compute new pair of (swUp, pc)
   if maxDiffSw > diffTolSw
      inx = maxInxSw;
   else
      inx = maxInxFf;
   end
   ffVal = mean(ffVec([inx, inx+1]));
   ff = ffVal*ones(G.cells.num,1);
   sw = fluid.fracFlowInv(ff);
   swVal = sum(sw.*pv) / pvTot;
   
   % Store value and keep sorted order
   ffVec = [ffVec(1:inx); ffVal; ffVec(inx+1:end-1)];
   swU   = [  swU(1:inx); swVal;   swU(inx+1:end-1)];
   
end

assert( (all(diff(ffVec) > 0) || all(diff(ffVec) < 0) ) && ...
   ( all(diff(swU) > 0) || all(diff(swU) < 0) ), ...
   'Fractional flow curve is not invertible');

% For a given upscaled water saturation, return fractional flow
ffFun = @(sWU) interp1(swU, ffVec, sWU);

end







