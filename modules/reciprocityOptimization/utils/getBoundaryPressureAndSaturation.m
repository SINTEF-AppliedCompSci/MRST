function newbc = getBoundaryPressureAndSaturation(model, states, bc)
%     p = model.AutoDiffBackend.convertToAD(zeros(numel(bc.face), 1), propsRes.pressure);
G = model.G;
T = model.operators.T;
if ~isempty(bc)
    is_pressure = reshape(strcmpi(bc.type, 'pressure'), [], 1);
    is_flux     = reshape(strcmpi(bc.type, 'flux'), [], 1);
    is_sat     = reshape(strcmpi(bc.type, 's'), [], 1);
    bdf = bc.face;
    cells = sum(model.G.faces.neighbors(bc.face,:), 2);
else
    bdf = boundaryFaces(G);
    cells = sum(model.G.faces.neighbors(bdf,:), 2);
    no_flux = bdf;
    is_pressure = [];
    is_flux = [];
end

for j = 1:length(states)
    state = states{j};
    % Get pressure at boundaries open to flow
    if any(is_pressure)
        if any(is_sat)
            newbc{j} = addBC([], bc.face(is_pressure), 'flux',        ...
                state.flux(is_pressure),'sat', bc.sat(is_pressure));
        else
            newbc{j} = addBC([], is_pressure, 'flux',        ...
                state.flux(is_pressure),'sat', state.s(cells(is_pressure)));
        end
    end
    if any(is_flux)
        % Flux BCs requires reconstruction. Assuming simple
        G = model.G;
        q = -bc.value(is_flux);
        if any(strcmpi(G.type,'topSurfaceGrid'))
            dzbc = model.gravity(3)*(G.cells.z(cells) - G.faces.z(bc.face));
        else
            g    = model.getGravityVector();
            dz   = G.faces.centroids(bc.face,:) - G.cells.centroids(cells,:);
            dzbc = -dz*g';
        end
        dzbc = dzbc(is_flux);
        [totMob, rhoMob] = deal(0);
        nph = model.getNumberOfPhases();
        for i = 1:nph
            mobi = state.FlowProps.Mobility{i}(is_flux);
            rhoi = state.PVTProps.Density{i}(is_flux);
            totMob = totMob + mobi;
            rhoMob = rhoMob + mobi.*rhoi;
        end
        fpress = (-q./T + rhoMob.*dzbc)./totMob;

        p  = xr.pressure;
   fp =  ...
          accumarray(G.cells.faces(:,1), (p(cellNo)+grav).*T, [G.faces.num,1])./ ...
          accumarray(G.cells.faces(:,1), T, [G.faces.num,1]);
        if any(is_sat)
            newbc{j} = addBC([], is_flux, 'pressure',        ...
                fpress,'sat', bc.sat(is_flux));
        else
            newbc{j} = addBC(newbc{j}, bdf, 'flux',        ...
            fpress,'sat', [state.s(cells(is_flux),:)]);
        end
    end
    if any(no_flux)
        if any(strcmpi(G.type,'topSurfaceGrid'))
            dzbc = model.gravity(3)*(G.cells.z(cells) - G.faces.z(bdf));
        else
            g    = model.getGravityVector();
            dz   = G.faces.centroids(bdf,:) - G.cells.centroids(cells,:);
            dzbc = -dz*g';
        end
        [totMob, rhoMob] = deal(0);
        nph = model.getNumberOfPhases();
        for i = 1:nph
            mobi = state.FlowProps.Mobility{i}(cells);
            rhoi = state.PVTProps.Density{i}(cells);
            totMob = totMob + mobi;
            rhoMob = rhoMob + mobi.*rhoi;
        end
        fpress = (rhoMob.*dzbc)./totMob;
        newbc{j} = addBC(newbc{j}, bdf, 'flux',        ...
            fpress,'sat', [state.s(cells,:)]);
    end
end