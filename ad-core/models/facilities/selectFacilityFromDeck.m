function facility = selectFacilityFromDeck(deck, model)
    if isempty(model.FacilityModel)
        facility = FacilityModel(model);
    else
        facility = model.FacilityModel;
    end
    % Vertical lift tables
    if isfield(deck.RUNSPEC, 'VFPIDIMS') || ...
       isfield(deck.SCHEDULE.control(1), 'VFPPROD') || ...
       isfield(deck.SCHEDULE.control(1), 'VFPINJ')
        facility.VFPTablesInjector = {};
        facility.VFPTablesProducer = {};
        for cNo = 1:numel(deck.SCHEDULE.control)
            vfpinj = deck.SCHEDULE.control(cNo).VFPINJ;
            for tno = 1:numel(vfpinj)
                if ~isempty(vfpinj{tno})
                    if numel(facility.VFPTablesInjector) >= tno
                        assert(isempty(facility.VFPTablesInjector{tno}), ...
                            'MRST does not presently support overwriting VFP tables');
                    end
                    facility.VFPTablesInjector{tno} = VFPTable(vfpinj{tno});
                end
            end
            
            vfpprod = deck.SCHEDULE.control(cNo).VFPPROD;
            for tno = 1:numel(vfpprod)
                if ~isempty(vfpprod{tno})
                    if numel(facility.VFPTablesProducer) >= tno
                        assert(isempty(facility.VFPTablesProducer{tno}), ...
                            'MRST does not presently support overwriting VFP tables');
                    end
                    facility.VFPTablesProducer{tno} = VFPTable(vfpprod{tno});
                end
            end

        end
    end
end