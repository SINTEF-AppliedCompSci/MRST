classdef CompositionalBrineFluid < CompositionalMixture
%CompositionalBrineFluid   Compositional mixture for H2O/NaCl (brine) systems
%
%   This class represents a compositional fluid mixture for geothermal
%   simulations involving water and salt (NaCl). It extends the
%   CompositionalMixture class to include molecular diffusivity and is
%   intended for use in geothermal and hydrothermal reservoir models where
%   brine composition and transport are important.
%
%   Properties:
%     molecularDiffusivity - Molecular diffusivity of salt in water (m^2/s)
%
%   Example usage:
%     fluid = CompositionalBrineFluid({'H2O', 'NaCl'}, [18.015, 58.44], 1e-9);
%
%   See also: CompositionalMixture, GeothermalModel

    properties
        molecularDiffusivity;
    end
    
    methods
        %-----------------------------------------------------------------%
        function fluid = CompositionalBrineFluid(names, molarMass, molecularDiffusivity, varargin)
        % Constructor for brine fluid mixture
        %
        %   fluid = CompositionalBrineFluid(names, molarMass, molecularDiffusivity, ...)
        %
        %   Parameters:
        %     names                - Cell array of component names (e.g., {'H2O', 'NaCl'})
        %     molarMass            - Vector of molar masses for each component
        %     molecularDiffusivity - Scalar or vector of molecular diffusivities (m^2/s)
        %     ...                  - Additional options passed to CompositionalMixture
        %
        %   Returns:
        %     fluid - Instance of CompositionalBrineFluid
        %
        %   This constructor sets up a compositional brine fluid for use in geothermal
        %   simulations, with support for molecular diffusion of salt.
            fluid = fluid@CompositionalMixture(names, nan, nan, nan, nan, molarMass, varargin{:});
            fluid.molecularDiffusivity = molecularDiffusivity;
        end
    end
    
end

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