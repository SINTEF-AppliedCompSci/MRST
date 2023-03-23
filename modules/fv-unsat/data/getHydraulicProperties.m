function soil = getHydraulicProperties(soil_name)
% Catalog of unsaturated hydraulic properties for different soils
%
% SOURCES: 
% [1] : https://doi.org/10.1029/WR024i005p00755
% [2] : https://doi.org/10.1029/WR026i007p01483
%
% SYNOPSIS:
%   soil = getHydraulicProperties(soil_name)
%
% PARAMETERS:
%   soil_name   - String, name of the soil. The catalog includes the
%                 following types of soil:
%       
%       Clay [1]                -> 'clay' 
%       Clay loam [1]           -> 'clayLoam'
%       Loam [1]                -> 'loam'
%       Loamy sand [1]          -> 'loamySand'
%       Silt [1]                -> 'silt'
%       Silt loam [1]           -> 'siltLoam'
%       Silt clay [1]           -> 'siltClay'
%       Silty clay loam [1]     -> 'siltyClayLoam'
%       Sand [1]                -> 'sand'
%       Sandy clay [1]          -> 'sandyClay'
%       Sandy clay loam [1]     -> 'sandyClayLoam'
%       Sandy loam [1]          -> 'sandyLoam'
%       New Mexico sample [2]   -> 'newMexSample'   
%
%  RETURNS:
%   soil        - Structure, containg the following properties:
%
%       theta_s [-]       -> water content at saturated conditions
%       theta_r [-]       -> residual water content
%       K_s [m/s]         -> saturated hydraulic conductivity     
%       alpha [1/m]       -> van Genuchten-Mualem equation parameter
%       n [-]             -> van Genuchten-Mualem equation parameter
%
%  EXAMPLE:
%   >> soil = getHydraulicProperties('sand')
% 
%   soil = 
% 
%       struct with fields:
% 
%       theta_s: 0.4300
%       theta_r: 0.0450
%           K_s: 8.2500e-05
%         alpha: 14.5000
%             n: 2.6800
%

%{
Copyright 2018-2020, University of Bergen.

This file is part of the fv-unsat module.

fv-unsat is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

fv-unsat is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this file.  If not, see <http://www.gnu.org/licenses/>.
%} 

% Soil structure
soil = struct('theta_s', [], 'theta_r', [], 'K_s', [], 'alpha', [], 'n', []);

% Catalog of soils
catalog = {'clay', 'clayLoam', 'loam', 'loamySand', 'silt', ...
           'siltLoam', 'siltyClay', 'siltyClayLoam', 'sand', ...
           'sandyClay', 'sandyClayLoam', 'sandyLoam', 'newMexSample'};
theta_s = [0.38, 0.41, 0.43, 0.41, 0.46, 0.45, 0.36, 0.43, 0.43, 0.38, 0.39, 0.41, 0.368];
theta_r = [0.068, 0.095, 0.078, 0.057, 0.034, 0.067, 0.070, 0.089, 0.045, 0.100, 0.100, 0.065, 0.102];
K_s     = [0.20, 0.26, 1.04, 14.59, 0.25, 0.45, 0.02, 0.07, 29.70, 0.12, 1.31, 4.42, 33.1920];
alpha   = [0.008, 0.019, 0.036, 0.124, 0.016, 0.020, 0.005, 0.010, 0.145, 0.027, 0.059, 0.075, 0.0335];
n       = [1.09, 1.31, 1.56, 2.28, 1.37, 1.41, 1.09, 1.23, 2.68, 1.23, 1.48, 1.89, 2];

% Obtain properties
switch soil_name 
    case catalog
        idx = find(strcmp(catalog, soil_name));
        soil.theta_s = theta_s(idx);
        soil.theta_r = theta_r(idx);
        soil.K_s = K_s(idx) * centi * meter / hour;
        soil.alpha = alpha(idx) / (centi * meter);
        soil.n = n(idx);
    otherwise
        error('Soil type not found. See documentation for available options.')
end       
