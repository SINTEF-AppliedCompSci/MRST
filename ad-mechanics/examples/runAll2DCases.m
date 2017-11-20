%%  Example of poelasticity simulation on 2D grids
%
% The implementation of the poroelasticty solver is done in such a way that it
% is easy to switch between different fluid models and different solving
% strategy (fully coupled versus splitting). In this example, we run all the
% combination of these choices for a 2D case. See run2DCase for a presentation
% of the test case.
%
% Options that are run are
%
% - fluid model : 
%
%   * 'water'     : single phase model
%   * 'oil water' : two phases model
%   * 'blackoil'  : Three phases, blackoil model
%
% - solver : 
%
%   * 'fully coupled'          : fully coupled solver
%   * 'fixed stress splitting' : solver using fixed stress splitting 
%


clear opt

opt.fluid_model        = 'water';
opt.cartDim            = [100, 10];
opt.L                  = [100, 10];
opt.method             = 'fully coupled';
opt.bc_case            = 'bottom fixed';

opt.nonlinearTolerance = 1e-6;
opt.splittingTolerance = 1e-6;
opt.verbose            = false;
opt.splittingVerbose   = false;

% write intro text for each case
writeIntroText = @(opt)(fprintf('\n*** New simulation\n* fluid model : %s\n* method : %s\n\n', opt.fluid_model, opt.method));

%%  water cases
opt.fluid_model = 'water';

opt.method      = 'fully coupled';
writeIntroText(opt);
[model, states] = run2DCase(opt);
plot2DPoroExample(model, states, opt);


opt.method      = 'fixed stress splitting';
writeIntroText(opt);
[model, states] = run2DCase(opt);
plot2DPoroExample(model, states, opt);

%%  Two phase oil water phases cases
opt.fluid_model = 'oil water';


opt.method      = 'fully coupled';
writeIntroText(opt);
[model, states] = run2DCase(opt);
plot2DPoroExample(model, states, opt);


opt.method      = 'fixed stress splitting';
writeIntroText(opt);
[model, states] = run2DCase(opt);
plot2DPoroExample(model, states, opt);

%%  Three phases Black-Oil phases cases
opt.fluid_model = 'blackoil';


opt.method      = 'fully coupled';
writeIntroText(opt);
[model, states] = run2DCase(opt);
plot2DPoroExample(model, states, opt);


opt.method      = 'fixed stress splitting';
writeIntroText(opt);
[model, states] = run2DCase(opt);
plot2DPoroExample(model, states, opt);
