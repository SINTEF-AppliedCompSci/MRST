function [eqs, bq, sol] = getWellContributionsBO(W, pBH, qs, p, rho, b, rs, m, sol)
%Undocumented Utility Function

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

   nPerf = cellfun(@numel, {W.cells})';
   perf2well = rldecode((1:numel(W))', nPerf);
   Rw = sparse((1:numel(perf2well))', perf2well, 1, numel(perf2well), numel(W));
   Tw = vertcat(W(:).WI);

   % compute drawdown
   Hw = cell(numel(W),1);
   for wnr = 1:numel(W)
      nperf = numel(W(wnr).cells);
      if ~isfield(W(wnr), 'topo')||isempty(W(wnr).topo)
         W(wnr).topo = [(0:(nperf-1))', (1:nperf)'];
      end
      if(~isfield(sol(wnr),'alpha'))
        sol(wnr).alpha=[1 0 0];
      end
      if(isempty(sol(wnr).alpha))  
        sol(wnr).alpha=[1 0 0];
      end
            
      
      Hw{wnr} = computeWellHead(W(wnr), sol(wnr).alpha, subDb(rho, perf2well==wnr));
   end
   Hw = vertcat(Hw{:});
   drawdown = -(Rw*pBH+Hw) + p;
   sdd = sign(drawdown);

   %producing perforations (injecting are changed later on)
   bq = cell(1,3);
   for ph = 1:3
      bq{ph} = -Tw.*b{ph}.*m{ph}.*drawdown;
   end

   %calculate average wellbore mixture, non-trivial only for wells with crossflow
   perfInjInx = (sdd<0);
   [compi, crossFlowFlag] = calcCompi(W, bq, qs, b, rs, perfInjInx, perf2well);

   %calculate mobilities for injecting connections
   mtInj = m{1}(perfInjInx)+m{2}(perfInjInx)+m{3}(perfInjInx);
   for ph = 1:3
      m{ph}(perfInjInx) = compi{ph}(perfInjInx).*mtInj;
   end

   %recalculate bq (could also just update injecting perforations ...)
   for ph = 1:3
      bq{ph}  = -Tw.*b{ph}.*m{ph}.*drawdown;
   end

   %well equations
   eqs{1} = -Rw'*bq{1} + qs{1};
   eqs{2} = -Rw'*bq{2} + qs{2};
   eqs{3} = -Rw'*(bq{3} + rs.*bq{2}) + qs{3};
   %check limits and update comtrols
   % [withinLims, W]
   % checkLims(W, pBH, qWs, qOs, qGs)
   %[sol,sol] = checkLims(sol, pBH, qs{:});
   [limW, W] = checkLims(W, pBH, qs{:});                            %#ok
   %eqs{4}   = handleBC(sol, pBH, qs{:});
   eqs{4}   = handleBC(W, pBH, qs{:});

   % finally explicit update of segment phase distribution and add crossflow
   % flag and connection pressure to sol
   alpha = calcAlpha(W, bq, b, rs, perf2well);

   % figure(1);
   % nn = find(perf2well == 1);
   % plot(nn, alpha{1}(:,1), '*-', nn, compi{1}(nn));
   % figure(2)
   % plot(nn, b{1}(nn));

   pBHDb = double(pBH);
   for wnr = 1:numel(W)
      sol(wnr).crossflow = crossFlowFlag(wnr);
      sol(wnr).cpressure = pBHDb(wnr) + Hw(perf2well==wnr);
      sol(wnr).alpha     = alpha{wnr};
   end
end

%--------------------------------------------------------------------------

function sc = subDb(c, inx)
   if iscell(c)
      if isa(c{1}, 'ADI')
         assert (all(cellfun(@(x) isa(x, 'ADI'), c)), ...
                 'All array elements must be ADI if first element is ADI');

         sc = cellfun(@(x)x.val(inx), c, 'UniformOutput', false);
      else
         sc = cellfun(@(x)x(inx), c, 'UniformOutput', false);
      end
   else
      sc = c.val(inx);
   end
end
