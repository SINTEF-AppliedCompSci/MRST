function bc = BoudaryConditions(model, scheduleRow)
%
% DESCRIPTION: define the boundary conditions for the simulation
%
% SYNOPSIS:
%   bc = BoudaryConditions(model,scheduleRow)
%
% PARAMETERS:
%   model - struct containing following fields:
%   - experiment: saturation functions used for forward modeling
%   - grid
%   scheduleRow: the schedule line for which the boundary condition is set
%
% RETURNS:
%   bc: struct with the boundary conditions with the following fields:
%   - left
%   - right
%
% ----------------------------------
% (c) 2020-2022
% Siroos Azizmohammadi
% Omidreza Amrollahinasab
% Montanuniversität Leoben, Austria
% Chair of Reservoir Engineering
% https://dpe.ac.at/
% ----------------------------------
%
%%
G           = model.grid.G;  
process     = model.experiment.process;
processType = process.type;    
processName = lower(process.name);

bc = [];
switch lower(processType)
    case 'cent' % centrifuge
        rhoW       = model.experiment.fluid.rhoW.value;
        rhoO       = model.experiment.fluid.rhoNW.value;
        coreLength = model.experiment.geometry.length.value;
        switch lower(processName)
            case 'imbibition'
                deltap = scheduleRow(3) * mean(rhoW) * coreLength;
            case 'drainage' 
                deltap = scheduleRow(3) * mean(rhoO) * coreLength;                
        end
        SwL = 0; SoL = 1;
        SwR = 1; SoR = 0;
        bcl.sat = [SwL SoL];
        bcl.bhp = model.experiment.schedule.pout.value;
        bc      = pside(bc, G, 'xmin', bcl.bhp, 'sat', bcl.sat);
        bc.left = bcl;
    
        % constant pressure boundary condition on the right, bcr
        bcr.sat  = [SwR SoR];
        bcr.bhp  = model.experiment.schedule.pout.value;
        bc       = pside(bc, G, 'xmax', bcr.bhp  + deltap, 'sat', bcr.sat);
        bc.right = bcr;
    otherwise
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