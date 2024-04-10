function setup = ensemble_base_problem_simple_gpsnet_model(varargin)
% Creates a GPSNET model with two-phase flow between two injectors and two
% producers.
%
% SYNOPSIS:
%   setup = ensemble_base_problem_simple_gpsnet_model('pn1', pv1, ...)
%   setup = ensemble_base_problem_simple_gpsnet_model(fullSetup, 'pn1', pv1, ...)
%
% DESCRIPTION:
%   
%   Very simple example that generates a 3D reservoir with a two-phase flow
%   problem. It has two injectors and two producers, of which one of the
%   injectors are not in a corner. May be used as a stand-alone example 
%   definition, or to construct an instance of `MRSTExample` as 
%   example = MRSTExample('ensemble_base_problem_3D_reservoir');
%   At the time of writing, the main purpose is to use this example for a
%   base in an example ensemble simulation.
%
% OPTIONAL PARAMETERS:
%   This example currently does not take any optional inputs
%
% RETURNS:
%   setup - test case with the following fields: name, description,
%      options, state0, model, schedule, and plotOptions.
%      If the optional input fullSetup (see synopsis) is false, the
%      returned setup only contains name, description, and options.
%
% SEE ALSO:
%   `TestCase`, `testcase_template`, `testSuiteTutorial`

%{
Copyright 2009-2024 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}


    error('Deprecated. Please use injector_to_producer_network from the network-model module instead');

 
end
