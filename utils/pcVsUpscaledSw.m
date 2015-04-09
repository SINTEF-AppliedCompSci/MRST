function [pcFun, swUMin, swUMax] = pcVsUpscaledSw(G, rock, fluid, ...
   varargin)
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

opt = struct(...
   'nPointsMin',   30, ...
   'nPointsMax',   50, ... %550, ... % Overrules difftols (breaks after nPointsMax)
   'diffTolSw',  0.01, ... %0.005
   'diffTolPc',  0.01 ... %0.005 % Relative to pc scale
   );
opt = merge_options(opt, varargin{:});

assert(isfield(fluid, 'pcOWInv'), ...
   'The fluid structure must have the field ''pcOWInv''. ');

% Options
nPointsMin = opt.nPointsMin; % initial pc data points
nPointsMax = opt.nPointsMax; % max number of extra data points
diffTolSw  = opt.diffTolSw; % diff tolerence for upscaled sw
diffTolPc  = opt.diffTolPc; % diff tolerence for capillary pressure

% Find min and max pc values
pcMin = min( min(fluid.pcOW(zeros(G.cells.num,1))), ...
             min(fluid.pcOW( ones(G.cells.num,1))) );
pcMax = max( max(fluid.pcOW(zeros(G.cells.num,1))), ...
             max(fluid.pcOW( ones(G.cells.num,1))) );

% Allocate vectors
pcVec = nan(nPointsMax, 1);
swU   = nan(nPointsMax, 1); % Upscaled water sat

% Compute the initial set of equally spaces pc points
pcVec(1:nPointsMin) = linspace(pcMin, pcMax, nPointsMin)';
pv    = rock.poro.*G.cells.volumes;
pvTot = sum(pv);
for i = 1:nPointsMin
   pc = pcVec(i)*ones(G.cells.num, 1);
   sw = fluid.pcOWInv(pc); % inverse
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
   [maxDiffPc, maxInxPc] = max(abs(diff(pcVec)));
   
   % Scale maxDiffPc with pc span
   maxDiffPc = maxDiffPc/(pcMax - pcMin);
   
   if maxDiffSw < diffTolSw && maxDiffPc < diffTolPc
      if i == nPointsExtra+1
         break; % ok. Reached tolerence in last iteration
      else
         % Remove nans at end of vectors and break out of loop.
         last  = length(pcVec) - (nPointsExtra-i+1);
         pcVec = pcVec(1:last);
         swU  = swU(1:last);
      end
      break;
   elseif i == nPointsExtra+1
      if mrstVerbose
          warning(['Max number of points reached, but tolerence not met.\n'...
             'maxDiffSw=%1.2e, maxDiffPc=%1.2e.'], maxDiffSw, maxDiffPc);
      end
      break;
   end
   
   % Compute new pair of (swUp, pc)
   if maxDiffSw > maxDiffPc
      inx = maxInxSw;
   else
      inx = maxInxPc;
   end
   pcVal = mean(pcVec([inx, inx+1]));
   pc = pcVal*ones(G.cells.num,1);
   sw = fluid.pcOWInv(pc);
   swVal = sum(sw.*pv) / pvTot;
   
   % Store value and keep sorted order
   pcVec = [pcVec(1:inx); pcVal; pcVec(inx+1:end-1)];
   swU   = [  swU(1:inx); swVal;   swU(inx+1:end-1)];
   
end

assert( (all(diff(pcVec) > 0) || all(diff(pcVec) < 0) ) && ...
   ( all(diff(swU) > 0) || all(diff(swU) < 0) ), ...
   'Capillary pressure curve is not invertible');

% For a given upscaled water saturation, return capillary pressure
pcFun = @(sWU) interp1(swU, pcVec, sWU);
% sWUFun = @(pc) interp1(pcVec, swU, pc);

end

