function static = SetupStaticParams(model)
%
% DESCRIPTION: initialize the static parameters of the model
%
% SYNOPSIS:
%   static = SetupStaticParams(model)
%
% PARAMETERS:
%   model - struct containing following fields:
%   - simulation: time stepping and grid cells information
%   - experiment: saturation functions used for forward modeling
%   - plot: plotting options
%
% RETURNS:
%   static - structure with the static parameters of the model with fields:
%   - gauge: for the primary pressure tabs
%   - gauge_Pmid: in case we have measurements both from ends of the core
%                 and pressure tabs
%   - dt: simulation time steps
%   - fig: visualization info
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
if(isfield(model.simulation,'gaugeOff'))
    gaugeOff = model.simulation.gaugeOff.value;
    nCells   = model.simulation.nCells.value;
    length   = model.experiment.geometry.length.value;        
    dx = length / nCells;
    gauge.left = round(gaugeOff / dx) + 1;
    gauge.right = nCells - gauge.left + 1;
else
    gauge.left  = 1;
    gauge.right = model.grid.G.cells.num;
end   
static.gauge = gauge;    
if isfield(model.experiment.observation, "pressure_mid")
    pressure_mid = model.experiment.observation.pressure_mid;
    assert(isfield(pressure_mid,"gaugeOff"),"Enter the gauge off for the mid pressure")
    gaugeoff_Pmid = pressure_mid.gaugeOff.value;
    nCells   = model.simulation.nCells.value;
    length   = model.experiment.geometry.length.value;        
    dx = length / nCells;
    gauge_Pmid.left = round(gaugeoff_Pmid / dx) + 1;
    gauge_Pmid.right = nCells - gauge_Pmid.left + 1;
    static.gauge_Pmid = gauge_Pmid;
end
static.dt  = model.simulation.timeStep.value;   
%---------------------------------------
fig.title    = 'Visualization';
if(isfield(model.plot,'style'))
    if(strcmp(model.plot.style.inputStyle,'normal')), fig.style = 'normal'; end
    if(strcmp(model.plot.style.inputStyle,'docked')), fig.style = 'docked'; end
end    
fig.tag      = fig.title;
fig.subtitle = {'Flooding Experiment',...
                'Injection Rates',...
                'Saturation Front',...
                'Pressure Differential',...
                'Average Water Saturation',...
                'Production'};                
static.fig = fig;
%---------------------------------------