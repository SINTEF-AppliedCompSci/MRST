function [simRes, schedule, controls, output] = lineSearchAgr(...
    simRes, G, S, W, rock, fluid, schedule, controls, grad, ...
    objectiveFunction, stepSize, figProps, opt)
% Run agressive line search based on given gradient.
% Handles box-constraints and linear equality - and (probably not)
% inequality-constaints
% based on iteratively applying the constrints to the gradient until
% convergence.
%
% SYNOPSIS:
%  [simRes, schedule, controls, data] = lineSearch(...
%    simRes, G, S, W, rock, fluid, schedule, controls, grad, objectiveFunction, ...
%    stepSize, figProps, opt)
%
% DISCRIPTION:
%  1) Project gradient according to constraints iteratively until tollerance ConstTol is met
%     or max number of iterations MaxConstIts is met (returnes failure).
%     Constraints are applied in the following order:
%        i) Box const.
%       ii) Lin. UnEq. const.
%      iii) Lin. Eq. const.
%     The resulting norm of the projected gradient is used as stopping criteria
%
%  2) Line search along projected gradient:
%        Uses three points [x1 x2 x3] = [0 .5 1]*stepSize
%        (a) If projected gradient pgrad(x2)=pgrad(x3), then on boundary,
%            done
%        (b) If obj(x1)<obj(x2)<obj(x3), set stepSize=2*stepSize goto (a)
%        (c) If obj(x1)>obj(x2), set stepSize = .5*stepSize goto (a)
%        (d) Find max on quad-curve through (x1, x2, x3) done
%
%
% PARAMETERS:
%   G, S, W, rock, fluid - usual structures
%   simRes               - ...
%   schedule, controls   - ...
%   grad                 - gradien as given by computeGradient
%   objectiveFunction    - handle to objective function
%
%   'pn'/pv - List of 'key'/value pairs defining optional parameters.  The
%             supported options are:
%          - MaxPoints           : maximum number of points before to
%                                  terminate line search algorithm
%          - LinSrchTol          : line search relative tolerance
%          - StepSize            : Step size
%          - ConstTol            : Relative contraint satisfaction tol
%          - MaxConstIts         : max number of its for constraint satisfaction
%          - VerboseLevel        : amount of output to screen
%
% RETURNS:
%   simRes     - Results for 'best' forward simulation
%   schedule   -
%   controls   -
%   data       - structure with fields
%                value         : objective function value for best run
%                relNormGrad   : relative norm of gradient |du|/|objective|
%                success       : whether or line search procedure succeeded
%                fraction      : optimal fraction obtained during the line
%                                search
%
%
% SEE ALSO:

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


verboseLevel = opt.verboseLevel;
plotting     = ~isempty(figProps);

if verboseLevel >= 0, fprintf('\n********* Starting line search ***********\n'); end

%--------------------------------------------------------------------------
% Initialize
x    = [0; 1/2; 1];
vals = -inf*ones(3, 1);

obj  = objectiveFunction(G, S, W, rock, fluid, simRes);
vals(1) = obj.val;
cp      = [0, vals(1)];
uInit  = [controls.well.values]';

output.value       = obj.val;
output.success     = true;
output.fraction    = 1;
output.stepSize    = stepSize;

continueSearch = true;

% ----------- Scale gradient wrt control-types  ---------------------------
[grad, scaleFacs] = scaleGradient(grad, simRes, controls, obj);
du = [grad{:}];

% -------- First determine if norm of projected gradient is less than tol -
projGrad = @(duCur)projectGradient(controls, duCur, ...
                                   'ConstTol', opt.constTol, ...
                                   'MaxConstIts', opt.maxConstIt);


pdu       = projGrad(du);
normPdu   = norm(pdu ./ scaleFacs, Inf);
output.normGrad = normPdu;
if normPdu < opt.gradTol
    fprintf('\nProjected gradient below tolerance. Done.');
    return;
end

if stepSize < 0    % guess stepSize
    stepSize = .5/max(max(abs(pdu./uInit)));

    if verboseLevel >= 0, fprintf('\nEstimated step size :  %8.6g', stepSize); end
end

% ----------- Create function for computing perturbed control -------------
compVal  = @(duCur)computeVal(simRes(1).resSol, G, S, W, rock, fluid, schedule, ...
                           controls, uInit, duCur, objectiveFunction,...
                           'ConstTol', opt.constTol, ...
                           'MaxConstIts', opt.maxConstIt, ...
                           'VerboseLevel', opt.verboseLevel);

% --- Check if pdu2==pdu3.
pointsOnBoundary = false;
du2     = du * x(2) * stepSize;
pdu2    = projGrad(du2);
du3     = du * x(3) * stepSize;
pdu3    = projGrad(du3);
if norm( pdu3-pdu2, Inf)/norm(pdu3, Inf) <= sqrt(eps)
    pointsOnBoundary = true;
end

% ---------- Compute points x2, x3 ----------------------------------------
vals(3) = compVal(du3);
cp = [cp; [x(3)*stepSize, vals(3)]]; plotif(figProps); plotif(figProps, cp);
if ~pointsOnBoundary
    vals(2) = compVal(du2);
    cp = [cp; [x(2)*stepSize, vals(2)]]; plotif(figProps, cp);
else
    vals(2) = vals(3);
    if vals(2)>vals(1)
        continueSearch = false;
    end
    cp = [cp; [x(2)*stepSize, vals(2)]]; plotif(figProps, cp, '-sg');
end

% ---------- Increase stepsize if val3>val2>val1 and recompute x3 -------------
if continueSearch && (vals(3)>vals(2)) && (vals(2)>vals(1))
    increaseStepSize = true;
    while increaseStepSize
        if verboseLevel >= 0, fprintf('\nIncreasing step size :  %8.6g', stepSize*2); end
        stepSize = stepSize * 2;
        vals(1:2) = vals(2:3);
        x(1)=1/4;
        du2 = du3; pdu2 = pdu3;
        du3    = du  * x(3) * stepSize;
        pdu3    = projGrad(du3);
        if norm( pdu3-pdu2, Inf)/norm(pdu3, Inf) <= sqrt(eps)
            pointsOnBoundary = true;
            continueSearch   = false;
            increaseStepSize = false;
            cp = [cp; [x(2)*stepSize, vals(2)]]; plotif(figProps, cp, '-sg');
        else
            vals(3) = compVal(du3);
            if vals(3)<vals(2)
                increaseStepSize = false;
                continueSearch   = false;
            end
            cp = [cp; [x(3)*stepSize, vals(3)]]; plotif(figProps, cp);
        end
    end
end

% ---------- Reduce stepsize if val1>val2 and recompute x2 -------------
if continueSearch && (vals(1)>vals(2))
    reduceStepSize = true;
    while reduceStepSize
        if verboseLevel >= 0, fprintf('\nReducing step size :  %8.6g', stepSize/2); end
        stepSize = stepSize/2;
        vals(3) = vals(2);
        du3 = du2; pdu3 = pdu2;
        du2    = du  * x(2) * stepSize;
        pdu2    = projGrad(du2);
        if norm( pdu3-pdu2, Inf)/norm(pdu3, Inf) < sqrt(eps)
            %do nothing, keep reducing
            cp = [cp; [x(2)*stepSize, vals(2)]]; plotif(figProps, cp, '-og');
        else
            vals(2) = compVal(du2);
            if vals(1) < vals(2)
                reduceStepSize  = false;
                continueSearch  = false;
            end
            cp = [cp; [x(2)*stepSize, vals(2)]]; plotif(figProps, cp);
        end
    end
end

% if not stuck on boundary find max on parabel
if ~pointsOnBoundary
    a = ([x.^2 x x.^0]\vals)/stepSize;
    xOpt  = -a(2)/(2*a(1));
    if (xOpt>0)&&(xOpt<1)
        duOpt = du  * xOpt * stepSize;
        valOpt = compVal(duOpt);
        cp = [cp; [xOpt*stepSize, valOpt]]; plotif(figProps, cp, '-om');
    else
        valOpt = -inf;
        cp = [cp; [xOpt*stepSize, valOpt]];
    end
else
    xOpt   = x(2);
    valOpt = vals(2);
    duOpt  = du2;
end

%--------------------------------------------------------------------------
if vals(2) > valOpt
    xOpt = x(2); duOpt = du2; valOpt = vals(2);
end
cp = [cp; [xOpt*stepSize, valOpt]]; plotif(figProps, cp, '-og');

pduOpt    = projGrad(duOpt);
uOpt      = uInit + pduOpt;
controls  = updateControls(controls, uOpt(:));
schedule  = updateSchedule(controls, schedule);
simRes    = runSchedule(simRes(1).resSol, G, S, W, rock, fluid, schedule, ...
                        'VerboseLevel', verboseLevel);

% Update data-struct
output.value       = valOpt;
output.success     = true;
output.fraction    = xOpt;
output.stepSize    = stepSize;
%output.relNormGrad = norm(pduOpt ./ scaleFacs, Inf);


function controls = updateControls(controls, u)
% update controls - strucure based on controll vector u
numControlWells = numel(controls.well);
numControls     = length(u);

U = reshape(u, numControlWells, numControls/numControlWells)';
for k = 1: numControlWells
    controls.well(k).values = U(:, k);
end

function [val, pdu] = computeVal(resSol, G, S, W, rock, fluid, schedule, ...
                                 controls, u, du, objectiveFunction, varargin)
opt     = struct('ConstTol',    1e-3, ...
                 'MaxConstIts',  100, ...
                 'VerboseLevel',   0);
opt     = merge_options(opt, varargin{:});

if opt.VerboseLevel >= 0, fprintf('\nForward simulation: ');tic; end
pdu       = projectGradient(controls, du, varargin{:});
uCur      = u + pdu;
controls  = updateControls(controls, uCur(:));
schedule  = updateSchedule(controls, schedule);
simRes    = runSchedule(resSol, G, S, W, rock, fluid, schedule, ...
                            'VerboseLevel', opt.VerboseLevel);
obj       = objectiveFunction(G, S, W, rock, fluid, simRes);
val       = obj.val;
 if opt.VerboseLevel >= 0
        tt = toc;
        fprintf('%4.2f sec. Objective value = %10.7e', tt, val);
 end

function [grad, scaleFacs] = scaleGradient(grad, simRes, controls, obj)
nSteps = numel(simRes) - 1;
cq = zeros(nSteps, 1); cp = zeros(nSteps, 1);
for step = 1 : nSteps
    ws  = simRes(step+1).wellSol;
    q   = cellfun(@sum, {ws.flux} );
    cq(step) = .5*(max(q) - min(q));
    p   = [ws.pressure];
    cp  = max(p) - min(p);
end
cq = mean(cq); cp = mean(cp);

controlTypes = {controls.well.type}';
isRate       = strcmp( controlTypes , 'rate');
fac          = ( isRate*cq + (~isRate)*cp );

nCSteps = nSteps;%controls.numControlSteps;
for cStep = 1 : nCSteps
    grad{cStep} = ( (fac.^2)/abs(obj.val) ) .* grad{cStep};
end

scaleFacs = repmat( fac, 1, nCSteps);

function [] = plotif(figProps, cp, lst)
if ~isempty(figProps)
    if nargin == 1
        cla(figProps.lsAxes); set(figProps.lsAxes, 'NextPlot', 'add');
    else
        scp = sortrows(cp,1);
        plot(figProps.lsAxes,scp(:,1), scp(:,2), '-*b');
        if nargin > 2
            plot(figProps.lsAxes, cp(end,1), cp(end,2), lst, 'LineWidth', 4);
        end
    end
    drawnow
end
