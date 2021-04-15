function [simRes, schedule, controls, output] = ...
                optimizeObjective(G, S, W, rock, fluid, resSolInit, ...
                schedule, controls, objectiveFunction, varargin)
% optimizeObjective -- Run whole optimization proccess using ad-hoc line
% search. A linkage with an external optimizer or use of matlabs optimizer
% toolbox (e.g., fmincon)is recommended.
%
% SYNOPSIS:
%  [simRes, schedule, controls, output] = optimizeObjective( ...
%                       G, S, W, rock, fluid, resSolInit, ...
%                       schedule, controls, objectiveFunction, varargin)
%
% PARAMETERS:
%   G, S, W, rock, fluid - usual structures
%   resSolInit           - initial 'solution' minimum containing field
%                          resSol.sw
%   schedule, controls   - ...
%   objectiveFunction    - handle to objective function
%
%   'pn'/pv - List of 'key'/value pairs defining optional parameters.  The
%             supported options are:
%   - gradTol       : stopping criterion for scaled gradien
%   - objChangeTol  : stopping criterion for realitive change in objective
%                     function
%   - stepSize      : initial step size
%   - VerboseLevel  : level of output to screen during progress (default: 0).
%                              < 0  : no output
%                                0  : minumal
%                                1  : quite a bit
%                                2  : loads
%   - PlotProgress  : whether or not to plot output during progress
%                            (default: true)
%   - OutputLevel   : level of output given in structure output (default: 0)
%                              < 0  = nothing (empty)
%                                0  = objective function value for every
%                                     iteration
%                                1  = ?
%
% RETURNS:
%   simRes     - Results for last (optimal) simulation
%   schedule   -
%   controls   -
%   outPut     - as described above
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


opt     = struct('gradTol',       1e-3, ...
                 'objChangeTol',  1e-4, ...
                 'constTol',     1e-10, ...
                 'stepSize',        -1, ...
                 'maxIt',           20, ...
                 'maxConstIt',     100, ...
                 'verboseLevel',     0, ...
                 'plotProgress',  true);
opt     = merge_options(opt, varargin{:});
stepSize     = opt.stepSize;
verboseLevel = opt.verboseLevel;
plotting     = opt.plotProgress;

%--------------------------------------------------------------------------
% Figure properties if plotting
figProps = initFig(plotting);
%--------------------------------------------------------------------------
objValues      = [];
iterationNum   = 0;
normGrad       = inf;
objChange      = inf;
%successLineSearch   = true;

 % Initial stepsize used in line search
schedules{1} = schedule;

while ( normGrad >= opt.gradTol) && (objChange >= opt.objChangeTol) && (iterationNum <= opt.maxIt)
    iterationNum = iterationNum + 1;
    if verboseLevel >= 0
        fprintf('\n********** STARTING ITERATION %3d ****************', ...
                iterationNum);
        fprintf('\nCurrent stepsize: %6.5f', stepSize);
    end

    % FORWARD SOLVE
    if iterationNum == 1   % For >1, forward sim is done in lineSearch
        if verboseLevel == 0, fprintf('\nForward solve %3d: ', iterationNum); tic; end;

        simRes = runSchedule(resSolInit, G, S, W, rock, fluid, schedule, ...
                             'VerboseLevel', verboseLevel);

        if verboseLevel == 0, tt = toc; fprintf('%6.3f sec. \n', tt);end

        obj       = objectiveFunction(G, S, W, rock, fluid, simRes);
        objValue  = obj.val;
        objValues = objValue;

        if verboseLevel >= 0, fprintf('Initial function value: %6.3f\n', objValue); end
        if plotting, plotProgress(figProps, G, simRes, controls, schedule, objValues); end
    end


    % ADJOINT SOLVE
    if verboseLevel == 0, fprintf('\nAdjoint solve %3d: ', iterationNum); tic; end;

    adjRes = runAdjoint(simRes, G, S, W, rock, fluid, schedule, controls, ...
                        objectiveFunction, 'VerboseLevel', verboseLevel);

    if verboseLevel == 0, tt = toc; fprintf('%6.3f sec. \n', tt);end


    % COMPUTE GRADIENT
    grad   = computeGradient(W, adjRes, schedule, controls);

    %Run agressive line-search
    [simRes, schedule, controls, lsData] = lineSearchAgr( ...
                     simRes, G, S, W, rock, fluid, schedule, controls, grad, ...
                     objectiveFunction, stepSize, figProps, opt);

    objValues  = [objValues; lsData.value];
    objChange  = 1 - abs(objValues(end-1)/objValues(end));
    normGrad   = lsData.normGrad;
    % update next stepsize (guess)
    stepSize          = 1.2*lsData.stepSize * lsData.fraction;

    schedules{iterationNum} = schedule;

    if verboseLevel >= 0
        fprintf('\n\n********** REPORT ITERATION %3d **********************', ...
                iterationNum);
        fprintf('\n*')
        fprintf('\n*  Obtained value               : %6.3e', objValues(end));
        fprintf('\n*  Change  (current/tol)        : %6.3e / %6.3e', objChange, opt.objChangeTol);
        fprintf('\n*  Gradient norm (current/tol)  : %6.3e / %6.3e', normGrad, opt.gradTol);
        fprintf('\n*')
        fprintf('\n********************************************************\n')
    end
    plotProgress(figProps, G, simRes, controls, schedule, objValues);
end
output.values = objValues;
output.schedules = schedules;

function [] = plotProgress(figProps, G, simRes, controls, schedule, vals)
if ~isempty(figProps)
    plot(figProps.objAxes, vals,'-*');

    if isfield(figProps, 'satAxes')
        axes(figProps.satAxes);cla;
        plotCellData(G, 1- simRes(end).resSol.s(:,1));
        axis tight;caxis([0 1]);
    end
    if isfield(figProps, 'rateAxes')
        types  = horzcat(schedule(:).types)';
        ratInx = cellfun(@(x)strcmp(x, 'rate'), types);
        vals   = horzcat(schedule(:).values)';
        rates = vals; bhps = vals;
        rates(~ratInx) = NaN; bhps(ratInx) = NaN;
        rInx = rldecode((1:numel(schedule))', 2);
        xax  = [schedule.timeInterval]';
        %pr = @(varargin)plot(varargin{:},'Color', 'r');
        %pb = @(varargin)plot(varargin{:},'Color', 'b');
        pr = @(varargin)plot(varargin{:});
        pb = @(varargin)plot(varargin{:});
        plotyy(figProps.rateAxes, xax/day, rates(rInx,:)*day, ...
                                  xax/day, bhps(rInx,:)/barsa, ...
                                  pb, pr);
    end
    drawnow
end

function fp = initFig(plotting)
if ~plotting
    fp = [];
else
    fh = figure;
    set(fh, 'Name', 'Optimization status', 'NumberTitle', 'off', ...
            'Position', [5 8 338 814]);
    nplots = 4;

    objAx = subplot(nplots,1,1);
    set(objAx, 'NextPlot', 'add')
    title('Objective function evolution');
    lsAx  = subplot(nplots,1,2);
    title('Current line search')

    fp.figure  = fh;
    fp.objAxes = objAx;
    fp.lsAxes  = lsAx;

    if nplots > 2
        fp.rateAxes = subplot(nplots,1,3);
        title('Rates/BHPs');
    end

    if nplots > 3
        fp.satAxes = subplot(nplots,1,4);
        title('End-time saturation');
    end


end











