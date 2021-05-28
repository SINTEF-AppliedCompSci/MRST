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
mrstModule add ad-mechanics

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
writeIntroText = @(opt)(fprintf('\n*** Start new simulation with\n* fluid model : %s\n* method : %s\n\n', opt.fluid_model, opt.method));

%%  water cases
opt.fluid_model = 'water';

opt.method      = 'fully coupled';
writeIntroText(opt);

optvals = cellfun(@(x) opt.(x), fieldnames(opt), 'uniformoutput', false);
optlist = reshape(vertcat(fieldnames(opt)', optvals'), [], 1);
[model, states] = run2DCase(optlist{:});
plot2DPoroExample(model, states, opt);


opt.method      = 'fixed stress splitting';
writeIntroText(opt);
optvals = cellfun(@(x) opt.(x), fieldnames(opt), 'uniformoutput', false);
optlist = reshape(vertcat(fieldnames(opt)', optvals'), [], 1);
[model, states] = run2DCase(optlist{:});
plot2DPoroExample(model, states, opt);

%%  Two phase oil water phases cases
opt.fluid_model = 'oil water';


opt.method      = 'fully coupled';
writeIntroText(opt);
optvals = cellfun(@(x) opt.(x), fieldnames(opt), 'uniformoutput', false);
optlist = reshape(vertcat(fieldnames(opt)', optvals'), [], 1);
[model, states] = run2DCase(optlist{:});
plot2DPoroExample(model, states, opt);


opt.method      = 'fixed stress splitting';
writeIntroText(opt);
optvals = cellfun(@(x) opt.(x), fieldnames(opt), 'uniformoutput', false);
optlist = reshape(vertcat(fieldnames(opt)', optvals'), [], 1);
[model, states] = run2DCase(optlist{:});
plot2DPoroExample(model, states, opt);

%%  Three phases Black-Oil phases cases
opt.fluid_model = 'blackoil';


opt.method      = 'fully coupled';
writeIntroText(opt);
optvals = cellfun(@(x) opt.(x), fieldnames(opt), 'uniformoutput', false);
optlist = reshape(vertcat(fieldnames(opt)', optvals'), [], 1);
[model, states] = run2DCase(optlist{:});
plot2DPoroExample(model, states, opt);


opt.method      = 'fixed stress splitting';
writeIntroText(opt);
optvals = cellfun(@(x) opt.(x), fieldnames(opt), 'uniformoutput', false);
optlist = reshape(vertcat(fieldnames(opt)', optvals'), [], 1);
[model, states] = run2DCase(optlist{:});
plot2DPoroExample(model, states, opt);

%%
% <html>
% <p><font size="-1">
% Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.
% </font></p>
% <p><font size="-1">
% This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).
% </font></p>
% <p><font size="-1">
% MRST is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% </font></p>
% <p><font size="-1">
% MRST is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% </font></p>
% <p><font size="-1">
% You should have received a copy of the GNU General Public License
% along with MRST.  If not, see
% <a href="http://www.gnu.org/licenses/">http://www.gnu.org/licenses</a>.
% </font></p>
% </html>
