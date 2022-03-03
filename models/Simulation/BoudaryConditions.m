function bc = BoudaryConditions(model,scheduleRow)
    G           = model.grid.G;  
    process     = model.experiment.process;
    processType = process.type;    
    processName = lower(process.name);
    
    bc = [];
    if(strcmpi(processType,'cent')) % centrifuge
        rhoW       = model.experiment.fluid.rhoW.value;
        rhoO       = model.experiment.fluid.rhoNW.value;
        coreLength = model.experiment.geometry.length.value;
        switch processName
            case 'imbibition'
                deltap = scheduleRow(3) * mean(rhoW) * coreLength;
            case 'drainage' 
                deltap = scheduleRow(3) * mean(rhoO) * coreLength;                
        end
        SwL = 0; SoL = 1;
        SwR = 1; SoR = 0;
        bcl.sat = [SwL SoL];
%         bcl.bhp = model.experiment.schedule.pout.value;
        bcl.bhp = scheduleRow(3) * mean(rhoO) * (model.experiment.schedule.centRad.value - coreLength/2);
        bc      = pside(bc, G, 'xmin', bcl.bhp, 'sat', bcl.sat);
        bc.left = bcl;

        % constant pressure boundary condition on the right, bcr
        bcr.sat  = [SwR SoR];
        bcr.bhp  = model.experiment.schedule.pout.value;
        bcr.bhp = scheduleRow(3) * mean(rhoO) * (model.experiment.schedule.centRad.value + coreLength/2);
%         bc       = pside(bc, G, 'xmax', bcr.bhp  + deltap, 'sat', bcr.sat);
        % discussion with Holger: the boundary pressures for the containaer
        % in the centrifuge pressures are better to represent the real
        % case as the commented lines for a drainage CENT but no diffrence
        % is observed in the results
        bc       = pside(bc, G, 'xmax', bcr.bhp, 'sat', bcr.sat);
        bc.right = bcr;
    else
        qinj   = scheduleRow(3); bcl.qinj = qinj;
        fo_inj = scheduleRow(4); bcl.fo_inj = fo_inj;
        fw_inj = scheduleRow(5); bcl.fw_inj = fw_inj;
        qo_inj = scheduleRow(6); bcl.qo_inj = qo_inj;
        qw_inj = scheduleRow(7); bcl.qw_inj = qw_inj;
        SwL = fw_inj; SoL = fo_inj;
        SwR = 1 - fw_inj; SoR = 1 - fo_inj;
        bcl.sat = [SwL SoL];
        bc      = fluxside(bc, G, 'xmin', bcl.qinj, 'sat', bcl.sat);
        bc.left = bcl;

        % constant pressure boundary condition on the right, bcr
        bcr.sat = [SwR SoR];
        bcr.bhp = model.experiment.schedule.pout.value;
        bc = pside(bc, G, 'xmax', bcr.bhp, 'sat', bcr.sat);
        bc.right = bcr;
    end    
end