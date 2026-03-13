function dynamic = SetupDynamicParams(model)
%
% DESCRIPTION: initialize the dynamic parameters which will be saved during
%              the simulation
%
% SYNOPSIS:
%   dynamic = SetupDynamicParams(model)
%
% PARAMETERS:
%   model - struct containing following fields:
%   - state: first state of the model
%
% RETURNS:
%   dynamic: struct with the dynamic parameters to be saved during the
%   simulation
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
params.pDiff       = zeros(1,1);
params.pDiff_mid   = zeros(1,1);
params.qinj        = zeros(1,2);
params.Qinj        = zeros(1,2); 
params.qprod       = zeros(1,2);
params.Qprod       = zeros(1,2); 
params.qp_net      = zeros(1,2);
params.Qp_net      = zeros(1,2); 
params.PVI         = zeros(1,1); 

params.periodStart = []; 
params.periodEnd   = []; 
params.periodInterval   = [];
params.scheduleSteps    = zeros(1,1);
params.cumScheduleSteps = zeros(1,1);
params.SwAvg   = mean(model.state.s(:,1));
params.Sw_min   = 1; % set 1 for visualization
params.Sw_max   = 0; % set 0 for visualization
params.plotNo  = 0;
params.counter = 1;
params.fromIdx = 1;
params.toIdx   = 1;
params.schedule = [];
params.updateAxes = {1,3};
params.selectedTimeAxes = {2,4,5,6};
dynamic.states = []; 
dynamic.params = params;
dynamic.schedulereport = [];    