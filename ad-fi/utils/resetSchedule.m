function deck = resetSchedule(deck, smry)
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

nm = smry.getNms('TIME'); nm = nm{1};
tms = smry.get(nm, 'TIME', :); tms = tms(:);

if isfield(deck.RUNSPEC, 'SI')
   tms = convertFrom(tms, day);
end

tmr = [0; cumsum( deck.SCHEDULE.step.val )];

isRepStep = false(numel(tms), 1);

tmrCount = 1;
for k = 1:numel(tms)
    if tms(k)>tmr(tmrCount)*(1+10*eps)
        isRepStep(k) = true;
        tmrCount = tmrCount +1;
    end
end

repStep = isRepStep(2:end);
deck.SCHEDULE.step.val = diff(tms);
deck.SCHEDULE.step.repStep  = repStep;
deck.SCHEDULE.step.control  = deck.SCHEDULE.step.control(cumsum(repStep));
end
