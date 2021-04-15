function [W] = updateWells(W, scheduleStep)
% Update wells based on schedule time step

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


numWells   = numel(W);
multiscale = isfield(W(1), 'CS');

for k = 1 : numWells
    type      = scheduleStep.types{k};
    val       = scheduleStep.values(k);
    W(k).type = type;
    W(k).val  = val;
    if strcmp(type, 'bhp'),
      W(k).S.RHS.f = - val(ones([W(k).S.sizeB(1), 1]));
      W(k).S.RHS.h = 0;                % never used
      if multiscale
          numCoarseCells = length( W(k).coarseCells );
          W(k).CS.RHS.f = -W(k).val(ones([numCoarseCells, 1]));
          W(k).CS.RHS.h = 0;      % never used
      end
   elseif strcmp(type, 'rate'),
      W(k).S.RHS.f = zeros([W(k).S.sizeB(1), 1]);
      W(k).S.RHS.h = - val;
      if multiscale
          numCoarseCells = length( W(k).coarseCells );
          W(k).CS.RHS.f = zeros([numCoarseCells, 1]);
          W(k).CS.RHS.h = -W(k).val;
      end
   end
end
