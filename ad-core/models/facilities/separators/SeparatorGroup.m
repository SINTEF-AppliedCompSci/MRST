classdef SeparatorGroup
    properties
        phaseDestination
        separators
        topologicalOrder
        surfaceSeparator
        mode = 'moles'
    end
    
    properties (Access = protected)
    end
    
    methods
        function g = SeparatorGroup(separator, pressure, T, phaseDestination, varargin)
            if nargin == 1
                [pressure, T, phaseDestination] = deal([]);
            end
            nsep = numel(pressure);
            g.separators = cell(nsep, 1);
            for i = 1:nsep
                sep = separator;
                sep.pressure = pressure(i);
                sep.T = T(i);
                g.separators{i} = sep;
            end
            if size(phaseDestination, 1) > 0
                toNext = find(phaseDestination(:, 1) == 0);
                % Zero for liquid indicates to the next stage
                phaseDestination(toNext, 1) = toNext+1;
                % Final always goes to tank
                phaseDestination(end, 1) = 0;
                % Finally -1 means also goes to tank
                phaseDestination(phaseDestination(:, 1) == -1, 1) = 0;
                g.phaseDestination = phaseDestination;
                g.topologicalOrder = getTopologicalSort(g);
            end
            g = merge_options(g, varargin{:});
            if isempty(g.surfaceSeparator)
                % Standard conditions
                g.surfaceSeparator = separator;
            end
        end
        
        function [surfaceRates, rhoS, phaseMoleFractions] = getSurfaceRates(sg, model, componentMassRates)
            switch lower(sg.mode)
                case 'mass'
                    streamMass = reshape(componentMassRates, 1, []);
                    [surfaceRates, ~, phaseMoleFractions, rhoS] = sg.computeSurfaceRates(model, streamMass);
                case 'moles'
                    mw = cellfun(@(x) x.molarMass, model.Components);
                    streamMole = reshape(componentMassRates, 1, []);
                    moleConvert = any(mw ~= 1);
                    hasWater = model.water;
                    if moleConvert
                        for i = 1:numel(componentMassRates)
                            streamMole{i} = streamMole{i}./mw(i);
                        end
                    end
                    if hasWater
                        waterMass = componentMassRates{end};
                        streamMole = streamMole(1:end-1);
                    end

                    [surfaceRates, ~, phaseMoleFractions, molarDensity] = sg.computeSurfaceRates(model, streamMole);
                    rhoS = molarDensity;
                    if moleConvert
                        for ph = 1:numel(rhoS)
                            w = 0;
                            xy = phaseMoleFractions{ph};
                            for c = 1:numel(xy)
                                w = w + xy{c}.*mw(c);
                            end
                            rhoS{ph} = rhoS{ph}.*w;
                        end
                    end
                    if hasWater
                        rhoW = model.fluid.rhoWS(1);
                        rhoS = [repmat(rhoW, numelValue(surfaceRates{1}), 1), rhoS];
                        surfaceRates = [{waterMass./rhoW}, surfaceRates];
                    end
                otherwise
                    error('Unknowm mode %s', sg.mode);
            end
        end
        
        function [phaseVolumeStream, phaseStream, phaseFractions, density] = computeSurfaceRates(sg, model, stream, varargin)
            n = numel(sg.separators);
            sep = sg.surfaceSeparator;
            if n == 0
                % We have no separator, just the surface conditions
                switch lower(sg.mode)
                    case 'mass'
                        [phaseStream, phaseFractions, density] = sep.separateComponentMassStream(model, stream);
                    case 'moles'
                        [phaseStream, phaseFractions, density] = sep.separateComponentMoleStream(model, stream);
                end
                phaseVolumeStream = phaseStream;
                for i = 1:numel(phaseVolumeStream)
                    phaseVolumeStream{i} = phaseVolumeStream{i}./density{i};
                end
            else
                % We have one or more separators, perform full traversal
                surfaceMoleRate = sg.computeSurfaceMoleRates(model, stream, varargin{:});
                nph = numel(surfaceMoleRate);
                
                phaseFractions = surfaceMoleRate;
                phaseStream = cell(1, nph);
                phaseVolumeStream = cell(1, nph);
                density = cell(1, nph);
                ncomp = numel(phaseFractions{1});
                % Flash the final stream and take the total volume
                for ph = 1:nph
                    total = 0;
                    for c = 1:ncomp
                        total = total + phaseFractions{ph}{c};
                    end
                    phaseStream{ph} = total;
                    for c = 1:ncomp
                        phaseFractions{ph}{c} = phaseFractions{ph}{c}./total;
                    end
                    switch lower(sg.mode)
                        case 'mass'
                            [stream, ~, localdensity] = sep.separateComponentMassStream(model, surfaceMoleRate{ph});
                        case 'moles'
                            [stream, ~, localdensity] = sep.separateComponentMoleStream(model, surfaceMoleRate{ph});
                    end
                    density{ph} = localdensity{ph};
                    tmp = 0;
                    for i = 1:nph
                        tmp = tmp + stream{i}./localdensity{i};
                    end
                    phaseVolumeStream{ph} = tmp;
                end
            end
        end

        
        function [surfaceStreams, separatorStreams] = computeSurfaceMoleRates(sg, model, stream_mole, start)
            if nargin < 4
                % Starting separator, defaults to one
                start = 1;
            end
            assert(iscell(stream_mole));
            ncomp = numel(stream_mole);
            order = sg.topologicalOrder;
            % Number of stages
            n = numel(sg.separators);
            % Phase molefractions out from each separator
            out_fractions = cell(n, 1);
            % Phase mole streams
            out_molestreams = cell(n, 1);
            % Phase molar densities
            out_densities = cell(n, 1);
            
            dest = sg.phaseDestination;
            nph = size(dest, 2);
            % Total stream into each separator
            in_streams = cell(n, 1);
            % Initialize
            in_streams{start} = stream_mole;
            stream = cell(1, ncomp);
            [stream{:}] = deal(0);
            surfaceStreams = cell(1, nph);
            [surfaceStreams{:}] = deal(stream);
            for index = find(order==start):n
                sepNo = order(index);
                % Get separator and in-stream (will be correct since we are
                % following the topological order
                sep = sg.separators{sepNo};
                s = in_streams{sepNo};
                if isempty(s)
                    % Not active
                    continue
                end
                [ms, frac, dens] = sep.separateComponentMoleStream(model, s);
                % Send to liquid/vapor separators
                for ph = 1:nph
                    target = dest(sepNo, ph);
                    val = cell(1, ncomp);
                    for c = 1:ncomp
                        val{c} = frac{ph}{c}.*ms{ph};
                    end
                    if target == 0
                        % Accumulate into surface streams
                        for c = 1:ncomp
                            surfaceStreams{ph}{c} = surfaceStreams{ph}{c} + val{c};
                        end
                    else
                        % Otherwise accumulate into moles for target stream
                        if isempty(in_streams{target})
                            in_streams{target} = val;
                        else
                            for c = 1:ncomp
                                in_streams{target}{c} = in_streams{target}{c} + val{c};
                            end
                        end
                    end
                end
                % Store outputs
                out_molestreams{sepNo} = ms;
                out_fractions{sepNo} = frac;
                out_densities{sepNo} = dens;
            end
            if nargout > 1
                separatorStreams = struct('molestreams',    {out_molestreams},...
                                          'densities',      {out_densities}, ...
                                          'molefractions',  {out_fractions});
            end
        end
        
        function sorting = getTopologicalSort(sg)
            N = sg.phaseDestination;
            n = numel(sg.separators);
            % Final node -> surface conditions
            N(N == 0) = n+1;
            start = (1:n)';
            A = sparse([start, start], N, 1, n+1, n+1);
            if exist('digraph', 'file')
                G = digraph(A);
                sorting = toposort(G)';
            else
                require matlab_bgl
                sorting = topological_order(A);
            end
        end
    end
end

%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

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
