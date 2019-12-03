classdef SurfactantAdsorption < StateFunction
    properties
    end

    methods
        function prop = SurfactantAdsorption(model, varargin)
            prop@StateFunction(model, varargin{:});
            prop = prop.dependsOn({'surfactant', 'surfactantmax'}, 'state'); % check mechanism
        end

        function ads = evaluateOnDomain(prop, model, state)
            [cs, csmax] = model.getProps(state, 'surfactant', 'surfactantmax');
            fluid = model.fluid;
            ads  = effads(cs, csmax, fluid);
        end
    end
end

function y = effads(cs, csmax, f)
   if f.adsInxSft == 2
      y = f.surfads(max(cs, csmax));
   else
      y = f.surfads(cs);
   end
end

%{
Copyright 2009-2019 SINTEF Digital, Mathematics & Cybernetics.

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
