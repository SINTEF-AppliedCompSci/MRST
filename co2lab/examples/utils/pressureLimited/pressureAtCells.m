function obj = pressureAtCells(model, states, schedule,cells,pstep, varargin)
% states.pressure is a cell array.
% schedule is only used for time steps.
% penalty is a scalar.
% plim is a cell array.

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

   opt.ComputePartials = false;
   %opt.cells=[];
   %opt.step=[];
   opt.tStep = [];
   opt = merge_options(opt, varargin{:});
   assert(max(cells)<model.G.cells.num);
   %num_timesteps = numel(schedule.step.val);
   tSteps = opt.tStep;
   if isempty(tSteps)
      numSteps = numel(states);
      tSteps = (1:numSteps)';
      %dts = schedule.step.val;
   else
      assert(numel(tSteps) == 1);
      numSteps = 1;
      %dts = schedule.step.val(opt.tStep);
   end
   
   obj = repmat({[]}, numSteps, 1);
   
   for step = 1:numSteps  
        state = states{tSteps(step)}; %@@ +1?      
        p = state.pressure;
        % keep track of amount over or amount under plim at each time step
        if opt.ComputePartials
            sG = state.s(:,2);   % place holders
            sGmax = state.sGmax; % place holders
            nW = numel(schedule.control(1).W);
            pBHP = zeros(nW, 1); % place holders
            qGs = pBHP;          % place holders
            qWs = pBHP;          % place holders
            [p, ~, ~, ~, ~, ~] = initVariablesADI(p, sG, sGmax, qWs, qGs, pBHP); 
        end         
      obj{step}=double(tSteps(step)==pstep)*p(cells);
   end
end
