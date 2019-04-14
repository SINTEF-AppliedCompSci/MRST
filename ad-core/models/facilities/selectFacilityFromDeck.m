function facility = selectFacilityFromDeck(deck, model)
    % Pick FacilityModel from input deck
    if isempty(model.FacilityModel)
        if isa(model, 'GenericBlackOilModel')
            facility = ExtendedFacilityModel(model);
        else
            facility = FacilityModel(model);
        end
    else
        facility = model.FacilityModel;
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
                    end
                    facility.VFPTablesInjector{tno} = VFPTable(vfpinj{tno});
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