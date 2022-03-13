function example = makeRandomTraining(example, rScale, bhpScale, shutin)
% Random perturbation of an existing simulation schedule
%
% SYNOPSIS:
%    problem = makeRandomTraining(example, rScale, bhpScale, shutin)
%
% DESCRIPTION:
%    The function takes a simulation problem as input and makes a random
%    perturbation of the well controls so that a subsequent simulation will
%    excite a larger variation in reservoir states and well responses. Each
%    sequence of constant well controls is subdivided into groups of four
%    steps. Each of these groups, except for the first, is given a random
%    perturbation around the existing well control. This means that if a
%    given control is valid for less than four consecutive steps, it will
%    not be perturbed
%
% INPUT PARAMETERS:
%    problem - Self-contained packed simulation problem which can be passed
%              to a number of routines, e.g:
%                 - simulatePackedProblem: Simulate and store output.
%                 - getPackedSimulatorOutput: Retrieve results.
%                 - simulatePackedProblemBackground: Simulate in seperate
%                   thread.
%              The function will change the simulation schedule found in
%              the 'schedule' field of the problem.
%
%    rScale -  Scaling of rate controls, either a scalar value or a
%              function handle. Scalar values are converted into a function
%              handle of the form
%                rscale = @(x,r) (1 + 2*rscale*(x-.5))*r
%              and then applied to a uniformly distributed random number x
%              and the original rate r. To replace this perturbation by
%              another custom made one, the user can instead pass a
%              function handle to an appropriate function.
%
%    bhpScale - Same as the rScale parameter, but applied to bhp controls
%
%    shutin  - Boolean variable. Specifies whether the perturbations
%              should include shut in of wells.
%
% RETURNS:
%    problem - a self-contained problem
%
% EXAMPLE:
%    % Load a standard Norne test case and perturb the rates by 25% and the
%    % bottom-hole pressure controls by 5%
%    trueEx  = TestCase('norne_simple_wo');
%    trainEx = makeRandomTraining(trueEx, 0.25, 0.05, false);

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

schedule.step  = example.schedule.step;
control        = example.schedule.control;
[ctrlNo,nstep] = rlencode(example.schedule.step.control,1);
inds           = [0; cumsum(nstep)];

if ~isa(rScale,'function_handle')
    rScale = @(x,y) (1 + 2*rScale*(x-0.5))*y;
end
if ~isa(bhpScale,'function_handle')
    bhpScale = @(x,y) (1 + 2*bhpScale*(x-0.5))*y;
end
ctrlInd = 0;
rng(1);
for m=1:numel(ctrlNo)
    ctrlVals = ceil(0.25*(1:nstep(m))') + ctrlInd;
    schedule.step.control(inds(m)+1:inds(m+1)) = ctrlVals;
    
    % Preserve the first well control
    ctrlInd = ctrlInd+1;
    schedule.control(ctrlInd) = control(ctrlNo(m));

    % Perturb the remaining ones
    for n=ctrlInd+1:max(ctrlVals)
        schedule.control(n) = control(ctrlNo(m));
        W = schedule.control(n).W;
        for i=1:numel(W)
            switch W(i).type
                case 'rate'
                    W(i).val = rScale(rand,W(i).val);
                    if ~isempty(W(i).lims) && ~isinf(W(i).lims.rate)
                        W(i).lims.rate = W(i).val;
                    end
                case 'bhp'
                    W(i).val = bhpScale(rand,W(i).val);
                    if ~isempty(W(i).lims) && ~isinf(W(i).lims.bhp)
                        W(i).lims.bhp = W(i).val;
                    end
            end
            if shutin && (W(i).sign<0) && (rand(1)>0.9)
                W(i).status = false;
            end
        end
        schedule.control(n).W = W;
    end
    ctrlInd = max(ctrlVals);
end
example.schedule = schedule;
example.name = [example.name '_rand'];
