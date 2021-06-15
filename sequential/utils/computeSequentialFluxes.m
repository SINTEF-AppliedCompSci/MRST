function [q_phase, q_components] = computeSequentialFluxes(state, G, vT, T, mob, rho, components, upstr, upwindType)
%Undocumented Utility Function

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

    nph = numel(mob);
    ncomp = numel(components);
    
    [flag_v, flag_g] = getSaturationUpwind(upwindType, state, G, vT, T, mob, upstr);
    isEqual = strcmpi(upwindType, 'potential') ||  all(flag_v(:) == flag_g(:));
    % Get components and mobility on the faces according to the first flag
    components_f = upstreamWeightCompositions(components, upstr, flag_v, nph, ncomp);
    mob_f = upstreamWeightPhases(mob, upstr, flag_v, nph);
    if isEqual
        q_phase = computePhaseVolumetricFluxes(vT, T, mob_f, G);
        if nargout > 1
            q_components = computeComponentFluxes(q_phase, components_f);
        end
    else
        % Different upwinding for potentials and viscous flux
        q_visc = computePhaseVolumetricFluxes(vT, T, mob_f);
        % Gravity part
        mob_g = upstreamWeightPhases(mob, upstr, flag_g, nph);
        q_grav = computePhaseVolumetricFluxes(0, T, mob_g, G);
        
        q_phase = q_visc;
        for i = 1:nph
            q_phase{i} = q_phase{i} + q_grav{i};
        end
        if nargout > 1
            if strcmpi(upwindType, 'hybrid_combined')
                % Take volumetric flux to be fixed
                q_components = computeComponentFluxes(q_phase, components_f);
            else
                % Take each term separately
                q_components_v = computeComponentFluxes(q_visc, components_f);
                q_components_g = computeComponentFluxes(q_grav, components_f);
                q_components = q_components_v;
                for i = 1:ncomp
                    q_components{i} = q_components{i} + q_components_g{i};
                end
            end
        end
    end
end

function mob_f = upstreamWeightPhases(value, upstr, flag, nph)
    mob_f = cell(nph, 1);
    for i = 1:nph
        v = value{i};
        if numelValue(v) == size(flag, 1)
            % Already upwinded
            mob_f{i} = v;
        else
            % Cell-value, must be upwinded
            mob_f{i} = upstr(flag(:, i), v);
        end
    end
end

function components_f = upstreamWeightCompositions(components, upstr, flag, nph, ncomp)
    components_f = cell(ncomp, 1);
    [components_f{:}] = deal(cell(nph, 1));
    for i = 1:nph
        for j = 1:ncomp
            c = components{j}{i};
            if ~isempty(c)
                components_f{j}{i} = upstr(flag(:, i), c);
            end
        end
    end
end

function q_c = computeComponentFluxes(q, components)
    nph = numel(q);
    nc = numel(components);
    q_c = cell(nc, 1);
    [q_c{:}] = deal(zeros(size(value(q{1}))));
    for c = 1:nc
        for p = 1:nph
            component_in_phase = components{c}{p};
            if ~isempty(component_in_phase)
                q_c{c} = q_c{c} + q{p}.*component_in_phase;
            end
        end
    end
end

function q = computePhaseVolumetricFluxes(vT, T, mob_f, G)
    if nargin < 4
        G = [];
    end
    nph = numel(mob_f);
    q = cell(nph, 1);
    mob_t = zeros(size(vT));
    for i = 1:nph
        mob_t = mob_t + mob_f{i};
    end
    
    for i = 1:nph
        f = mob_f{i}./mob_t;
        
        if isempty(G)
            q{i} = f.*vT;
        else
            g = zeros(size(vT));
            for j = 1:nph
                if i == j
                    continue
                end
                g = g + mob_f{j}.*(G{i} - G{j});
            end
            q{i} = f.*(vT + T.*g);
        end
        
    end
end
