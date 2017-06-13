
%% Example: Coupling between fluid flow (blackoil model) and rock mechanics (elasticity)
% 
%

mrstModule add ad-mechanics ad-core ad-props ad-blackoil vemmech deckformat mrst-gui

%% Summary of the options

% opt.norne_case = 'mini Norne'; % 'full' or 'mini Norne'
%
% 'full'       : 7392 cells
% 'mini Norne' :  605 cells
%
% opt.bc_case    = 'bottom fixed'; % 'no displacement' or 'bottom fixed'
%
% 'no displacement' : All nodes belonging to external faces have displacement
%                     equal to zero
% 'bottom fixed'    : The nodes that belong to the bottom have zero
%                     displacement, while a given pressure is imposed on
%                     the external faces that are not bottom faces.
%
% opt.method     = 'fully coupled'; % 'fully coupled' 'fixed stress splitting'
%
% 'fully coupled'          : The mechanical and flow equations are solved fully couplde.
%
% 'fixed stress splitting' : The mechanical and flow equations are solved
%                            sequentially using a fixed stress splitting
%
% opt.fluid_model = 'oil water'; % 'blackoil' 'water' 'oil water'
%
% 'blackoil'     : blackoil model is used for the fluid (gas is injected, see
%                  schedule below)
% 'oil water'    : Two phase oil-water
% 'water' : water model is used for the fluid


opt.verbose          = true;
opt.splittingVerbose = true;
opt.norne_case       = 'mini Norne';
opt.bc_case          = 'bottom fixed';


setTitle = @(opt)(sprintf('%s, %s', opt.fluid_model, opt.method));
clear states;
i = 1;

%%  water cases

opt.fluid_model = 'water';

opt.method      = 'fully coupled';
[model, initState, schedule] = setupNorneExamples(opt);
[wsols{i}, states{i}] = simulateScheduleAD(initState, model, schedule);
figure
plotToolbar(model.G, states{i}, 'outline', true);
title(setTitle(opt));
view([7, 40]);
i = i + 1;

opt.method      = 'fixed stress splitting';
[model, initState, schedule] = setupNorneExamples(opt);
[wsols{i}, states{i}] = simulateScheduleAD(initState, model, schedule);
figure
plotToolbar(model.G, states{i}, 'outline', true);
title(setTitle(opt));
view([7, 40]);
i = i + 1;


%%  Two phase oil water phases cases

opt.fluid_model = 'oil water';

opt.method      = 'fully coupled';
[model, initState, schedule] = setupNorneExamples(opt);
[wsols{i}, states{i}] = simulateScheduleAD(initState, model, schedule);
figure
plotToolbar(model.G, states{i}, 'outline', true);
view([7, 40]);
i = i + 1;

opt.method      = 'fixed stress splitting';
[model, initState, schedule] = setupNorneExamples(opt);
[wsols{i}, states{i}] = simulateScheduleAD(initState, model, schedule);
figure
plotToolbar(model.G, states{i}, 'outline', true);
title(setTitle(opt));
view([7, 40]);
i = i + 1;


%%  Three phases Black-Oil phases cases

opt.fluid_model = 'blackoil';

opt.method      = 'fully coupled';
[model, initState, schedule] = setupNorneExamples(opt);
[wsols{i}, states{i}] = simulateScheduleAD(initState, model, schedule);
figure
plotToolbar(model.G, states{i}, 'outline', true);
title(setTitle(opt));
view([7, 40]);
i = i + 1;

opt.method      = 'fixed stress splitting';
opt.splittingTolerance = 1e-3;
[model, initState, schedule] = setupNorneExamples(opt);
[wsols{i}, states{i}] = simulateScheduleAD(initState, model, schedule);
figure
plotToolbar(model.G, states{i}, 'outline', true);
title(setTitle(opt));
view([7, 40]);
