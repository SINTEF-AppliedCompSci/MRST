function model = CreateFluid(model)
%
% DESCRIPTION: adds fluid model to the model struct
%
% SYNOPSIS:
%   model = CreateFluid(model)
%
% PARAMETERS:
%   model - struct containing following fields:
%   - experiment: saturation functions used for forward modeling
%   - simulation: time stepping and grid cells information
%
% RETURNS:
%   model - adds the the following fields to the struct:
%   - fluid: struct with the fluid model information
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
muW  = model.experiment.fluid.muW.value;
muO  = model.experiment.fluid.muNW.value;
rhoW = model.experiment.fluid.rhoW.value;
rhoO = model.experiment.fluid.rhoNW.value;
simulation = model.simulation;

fluid = initSimpleADIFluid('phases', 'WO', ...
                           'mu'    , [muW, muO]);
fluid.isIncomp = true;
fluid.rhoWS = rhoW;
fluid.rhoOS = rhoO;  
if (isfield(simulation,'bCells'))   
    krw_1 = @(sw) sw;
    krw_2 = @(sw) interpTable(model.satfun.sw_kr,model.satfun.krw,sw);
    kro_1 = @(so) so;
    kro_2 = @(so) interpTable(model.satfun.sw_kr,model.satfun.kro,1-so);
    fluid.krW = {krw_1, krw_2};
    fluid.krO = {kro_1, kro_2};
    pc_1 = @(s) 0;
    pc_2 = @(sw) interpTable(model.satfun.sw_pc,model.satfun.pc,sw);
    fluid.pcOW = {pc_1, pc_2};
else
    fluid.krW = @(sw) interpTable(model.satfun.sw_kr,model.satfun.krw,sw);
    fluid.krO = @(so) interpTable(model.satfun.sw_kr,model.satfun.kro,1-so);
    fluid.pcOW = @(sw) interpTable(model.satfun.sw_pc,model.satfun.pc,sw);      
end    
model.fluid = fluid;
