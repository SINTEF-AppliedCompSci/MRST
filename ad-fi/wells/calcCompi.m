function [compi, crossFlowFlag] = calcCompi(W, bq, qs, b, rs, perfInjInx, perf2well)
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

   crossFlowFlag = false(numel(W), 1);
   if isa(rs, 'ADI')
      cph   = double2ADI(zeros(numel(perf2well),1), rs);
   else
      cph = zeros(numel(perf2well), 1);
   end
   compi = {cph, cph, cph};
   for wnr = 1:numel(W)
      nperf = numel(W(wnr).cells);
      pinx  = find((perf2well == wnr));
      isInj = perfInjInx(pinx);
      if sum(isInj)==0||sum(~isInj)==0 %no cross-flow
         for ph=1:3
            compi{ph}(pinx) = ones(nperf,1)*W(wnr).compi(ph);
         end
      else %crossflow
         crossFlowFlag(wnr) = true;
         bqw = {bq{1}(pinx), bq{2}(pinx), bq{3}(pinx)};
         rsw = rs(pinx);
         qsw = {qs{1}(wnr), qs{2}(wnr), qs{3}(wnr)};
         bw  = {b{1}(pinx), b{2}(pinx), b{3}(pinx)};

         isPrd = ~isInj; %producing connections
                         %first compute inflow at surface conds
         inflows = { -isPrd'*bqw{1}, -isPrd'*bqw{2}, ...
                     -isPrd'*(bqw{3} + rsw.*bqw{2}) };
         for ph = 1:3
            inflows{ph} = inflows{ph} + qsw{ph}*(qsw{ph}>=0)*(W(wnr).compi(ph)>0);
         end

         %calculate inflow at perforation conditions
         qInj  = { inflows{1}./bw{1}, ...
                   inflows{2}./bw{2}, ...
                   (inflows{3} - inflows{2}.*rsw)./bw{3} };
         qInjt = qInj{1} + qInj{2} + qInj{3};

         %finally set compi
         negligeable = abs(double(qInjt))<1e-10/day;
         is_ADI = isa(qInjt, 'ADI');
         for ph=1:3
            if is_ADI
               compi_well = double2ADI(qInjt.val, qInjt); % dummy values
            else
               compi_well = zeros(numel(qInjt), 1);
            end
            if any(~negligeable)
               compi_well(~negligeable) = qInj{ph}(~negligeable)./qInjt(~negligeable);
            end
            if any(negligeable)
               compi_well(negligeable) = ones(numel(negligeable), 1)*W(wnr)*compi(ph);
            end
            compi{ph}(pinx) = compi_well;
         end

      end
   end

end



