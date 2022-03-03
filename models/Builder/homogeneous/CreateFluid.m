function model = CreateFluid(model)
% <keywords>
%
% Purpose : assign the fluid model interpolators to the model
%
% Syntax :
%
% Input Parameters :
%
% Return Parameters :
%
% Description :
%
% Author : 
%    Siroos Azizmohammadi
%    Omidreza Amrollahinasab
%
% History :
% \change{1.0}{09-Nov-2021}{Original}
%
% --------------------------------------------------
% (c) 2021, Siroos Azizmohammadi,
% Omidreza Amrollahinasab
% Chair of Reservoir Engineering, University of Leoben, Austria
% email: info@dpe.ac.at
% url: dpe.ac.at
% --------------------------------------------------
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
end
