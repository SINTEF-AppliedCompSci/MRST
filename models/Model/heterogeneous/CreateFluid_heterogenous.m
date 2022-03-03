function model = CreateFluid_heterogenous(model, f)
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
        fluid.krW = {};
        fluid.krO = {};
        fluid.pcOW = {};
        for i = 1: length(f)
            fluid.pcOW{1, i} = @(sw) interpTable(model.satfun.sw_pc./f(i),model.satfun.pc,sw);
            fluid.krW{1, i} = @(sw) interpTable(model.satfun.sw_kr,model.satfun.krw,sw);
            fluid.krO{1, i} = @(so) interpTable(model.satfun.sw_kr,model.satfun.kro,1-so);
        end
        fluid.pcOW{1, 1} = @(s) 0; fluid.pcOW{1, end} = @(s) 0;
        fluid.krW{1, 1} = @(sw) sw; fluid.krW{1, end} = @(sw) sw;
        fluid.krO{1, 1} = @(so) so; fluid.krO{1, end} = @(so) so;
    else
        fluid.krW = {};
        fluid.krO = {};
        fluid.pcOW = {};
        f(1) = []; f(end) = [];
        for i = 1: length(f)
            fluid.pcOW{1, i} = @(sw) interpTable(model.satfun.sw_pc./f(i),model.satfun.pc,sw);
            fluid.krW{1, i} = @(sw) interpTable(model.satfun.sw_kr,model.satfun.krw,sw);
            fluid.krO{1, i} = @(so) interpTable(model.satfun.sw_kr,model.satfun.kro,1-so);
        end
    end    
    model.fluid = fluid;
end
