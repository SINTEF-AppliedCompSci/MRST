function T = generateCoolPropTable(prop, propX, propY, rangeX, rangeY, varargin)
%Generate property table using CoolProp (see [1] and coolprop.org)
%
%   [1] Bell et. al., Pure and Pseudo-pure Fluid Thermophysical Property
%       Evaluation and the Open-Source Thermophysical Property Library
%       CoolProp. DOI: 10.1021/ie4033999

    % Optional arguments
    opt = struct('fluid'       , 'IF97::Water', ...
                 'numX'        , 500          , ...
                 'numY'        , 500          , ...
                 'save'        , false        , ...
                 'readFromDisk', true         );
    opt = merge_options(opt, varargin{:});
    % Adjust input
    if diff(rangeX) == 0, opt.numX = 1; end
    if diff(rangeY) == 0, opt.numY = 1; end
    % Get file and table names
    [fileName, tableName] = getTableNames(prop, propX, propY, rangeX, rangeY, opt.numX, opt.numY);
    % Check if table already exists
    fn = fullfile(mrstPath('geothermal'), 'eos', 'tables', fileName);
    if opt.readFromDisk && exist(fn , 'file')
        data = load(fn); T = data.T; return
    end
    % Make table
    [z, x, y] = makeTable(opt.fluid, prop, propX, propY, rangeX, rangeY, opt.numX, opt.numY);
    T = struct('name', tableName, ...
               'x'   , x        , ...
               'y'   , y        , ...
               'data', z        );
    if opt.save
        % Save results
        save(fullfile(mrstPath('geothermal'), 'eos', 'tables', fileName), 'T');
    end
end

%-------------------------------------------------------------------------%
function name = translateToCoolProp(name)
    switch lower(name)
        case 'pressure'
            name = 'P';
        case 'enthalpy'
            name = 'H';
        case 'internalenergy'
            name = 'U';
        case 'temperature'
            name = 'T';
        case 'viscosity'
            name = 'V';
        case 'density'
            name = 'D';
        case 'quality'
            name = 'Q';
        case 'conductivity'
            name = 'L';
        case 'heatcapacity'
            name = 'U';
        otherwise
            error('Unknown property name')
    end
end

%-------------------------------------------------------------------------%
function [z, x, y] = makeTable(fluid, prop, propX, propY, rangeX, rangeY, numX, numY)

    x = linspace(rangeX(1), rangeX(2), numX)';
    y = linspace(rangeY(1), rangeY(2), numY)';
    
    z = zeros(numX, numY);
    prop  = translateToCoolProp(prop);
    propX = translateToCoolProp(propX);
    propY = translateToCoolProp(propY);
    
    getProp = @(x, y) py.CoolProp.CoolProp.PropsSI(prop, propX, x, propY, y, fluid);
    for i = 1:numX
        for j = 1:numY
            try
                z(i,j) = getProp(x(i), y(j));
            catch me
                warning(me.identifier                                , ...
                        '%s. Interpolating to nearest valid neighbor', ...
                        me.message                                   );
                if ~isnan(z(max(i-1,1), j))
                    z(i,j) = z(max(i-1,1),j);
                elseif ~isnan(z(i, max(j-1,1)))
                    z(i,j) = z(i, max(j-1,1));
                else
                    z(i,j) = nan;
                end
            end
        end
    end    
end

%-------------------------------------------------------------------------%
function [fileName, tableName]= getTableNames(prop, propA, propB, rangeA, rangeB, numA, numB)
    tableName = sprintf('CoolProp table %s(%s[%d], %s[%d])', prop, propA, numA, propB, numB);
    fileName  = sprintf('%s_%s-%d-%d-%d_%s-%d-%d-%d.mat', prop, propA, numA, rangeA(1), rangeA(2), propB, numB, rangeB(1), rangeB(2));
end

%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

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