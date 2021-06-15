function states_seq = correctSequentialBHP(modelseq, states_seq, schedule)
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

    modelseq.transportModel.extraStateOutput = true;
    for i = 1:numel(states_seq)
        i
        W = schedule.control(schedule.step.control(i)).W;
        seq = states_seq{i};
        [~, seq] = transportEquationsNaturalVars(seq, seq, modelseq.transportModel, 1, schedule.control, 'resOnly', true);
        for j = 1:numel(W)
            w = W(j);
            q = sum(seq.wellSol(j).flux, 2);
            c = w.cells;
            WI = w.WI;
            mob = sum(seq.mob(c, :), 2);
            p = seq.pressure(c);

            if ~strcmpi(w.type, 'bhp')
                bhp_new = q./(WI.*mob) - seq.wellSol(j).cdp + p;
                bhp_new = sum(bhp_new.*abs(q))./sum(abs(q));
                states_seq{i}.wellSol(j).bhp = bhp_new;
            end
        end
    end
end
