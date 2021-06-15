function [sol, withinLims] = updateSwitchedWellControls(wellmodel, model, sol, pBH, q_s)
%Check for violated well limits and switch controls. 
%
% SYNOPSIS:
%   [sol, withinLims] = ...
%              updateSwitchedWellControls(wellmodel, model, sol, pBH, q_s)
%
% PARAMETERS:
%   wellmodel   - Simulation well model.
%   model       - Simulation model.
%   sol         - List of current well solution structures
%   pBH         - Vector of well bhps
%   q_s         - List of vectors of well component volume-rates 
%                 (surface conds) 
%
% RETURNS:
%   sol         - Well solution structures with updated fields 'type' and
%                 'val' in case of control switching
%   withinLims  - Logical vector with 'false' corresponding to wells where 
%                 limits were violated and switching has been perfomed.Ben
%
% SEE ALSO:
%   WellModel, setupWellControlEquation

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
W = wellmodel.W;
if isempty(W)||isempty(W(1).lims)
    withinLims = true(numel(sol),1);
else
    allowWellSignChange = wellmodel.allowWellSignChange;

    nwells     = numel(sol);
    withinLims = true(nwells,1);

    pBH   = double(pBH);
    q_s   = cell2mat( cellfun(@double, q_s, 'UniformOutput', false) );

    for wnr = 1:numel(sol)
        if isfield(W(wnr), 'status') && ~W(wnr).status
            % Inactive well, skip any limit checks
            continue
        end
        lims = W(wnr).lims;

        pBHw  = pBH(wnr);
        q_sw  = q_s(wnr,:);
        qt_sw = sum(q_sw);
        if ~isnumeric(W(wnr).lims)
            if ~allowWellSignChange
                lims.vrat = -inf;
            else
                lims.vrat = -inf;
            end
            if sol(wnr).sign > 0   % injector
                modes   = {'bhp', 'rate', 'rate'};
                flags = [pBHw > lims.bhp, qt_sw > lims.rate, qt_sw < lims.vrat];
            else            % producer
                modes   = {'bhp', 'orat', 'lrat', 'grat', 'wrat', 'vrat'};
                
                % insert dummy limits for missing fields
                missing_fields = {modes{~cellfun(@(x) isfield(lims, x), modes)}};
                for f = missing_fields
                   lims = setfield(lims, f{:}, -inf);
                end
                
                flags = [pBHw       < lims.bhp,  ...
                    q_sw(2)         < lims.orat, ...
                    q_sw(1)+q_sw(2) < lims.lrat, ...
                    q_sw(3)         < lims.grat, ...
                    q_sw(1)         < lims.wrat, ...
                    qt_sw           > -lims.vrat];
            end
        else
            modes = {};
            flags = false(numel(sol), 1);
            assert(isinf(lims))
        end
        %limits we need to check (all others than w.type):
        chkInx = ~strcmp(sol(wnr).type, modes);
        vltInx = find(flags(chkInx), 1);
        if ~isempty(vltInx)
            withinLims(wnr) = false;
            modes  = modes(chkInx);
            switchMode = modes{vltInx};
            fprintf('Well %s: Control mode changed from %s to %s.\n', sol(wnr).name, sol(wnr).type, switchMode);
            sol(wnr).type = switchMode;
            sol(wnr).val  = lims.(switchMode);
        end
    end
end
end


