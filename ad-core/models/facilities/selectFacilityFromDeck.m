function facility = selectFacilityFromDeck(deck, model)
% Pick FacilityModel from input deck

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

    if isempty(model.FacilityModel)
        if isa(model, 'GenericBlackOilModel')
            facility = GenericFacilityModel(model);
        else
            facility = FacilityModel(model);
        end
    else
        facility = model.FacilityModel;
    end
    if numel(deck.SCHEDULE.control) == 0
        return;
    end
    % Vertical lift tables
    if isfield(deck.RUNSPEC, 'VFPIDIMS') || ...
       isfield(deck.SCHEDULE.control(1), 'VFPPROD') || ...
       isfield(deck.SCHEDULE.control(1), 'VFPINJ')
        facility.VFPTablesInjector = {};
        facility.VFPTablesProducer = {};
        hasWarned = false;
        for cNo = 1:numel(deck.SCHEDULE.control)
            vfpinj = deck.SCHEDULE.control(cNo).VFPINJ;
            for tno = 1:numel(vfpinj)
                if ~isempty(vfpinj{tno})
                    if numel(facility.VFPTablesInjector) >= tno
                        if ~isempty(facility.VFPTablesInjector{tno}) && ~hasWarned
                            warning('MRST does not presently support overwriting VFP tables');
                            hasWarned = true;
                        end
                        facility.VFPTablesInjector{tno} = VFPTable(vfpinj{tno});
                    end
                end
            end
            vfpprod = deck.SCHEDULE.control(cNo).VFPPROD;
            for tno = 1:numel(vfpprod)
                if ~isempty(vfpprod{tno})
                    if numel(facility.VFPTablesProducer) >= tno
                        if ~isempty(facility.VFPTablesProducer{tno}) && ~hasWarned
                            warning('MRST does not presently support overwriting VFP tables');
                            hasWarned = true;
                        end
                    end
                    facility.VFPTablesProducer{tno} = VFPTable(vfpprod{tno});
                end
            end
        end
    end
end
