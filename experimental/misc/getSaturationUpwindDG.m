function [flagV, flagG, upCellsV, upCellsG, s_v, s_G] = getSaturationUpwindDG(disc, faces, x, T, flux, state, g, mob, sdof, rdof, xdof)
    % Explicit calculation of upstream cells (Bernier & Jaffre)
    % for each quadrature point x on each face in faces.
    %
    % PARAMETERS:
    %   faces - Faces we want upstream cells for
    %   x     - Upstram cells are computed for face(ix) at x(ix,:)
    %   sdof  - Saturation degrees of freedom
    %   state - For dofPos
    %   T, vT, g, mob - Inputs to getSaturationUpwind
    %
    % RETURNS:
    %   flag_v, flag_G       - Upstream flags for viscous and gravity
    %   upCells_v, upCells_G - Upstream cells for viscous and gravity
    %   s_v, sG              - Corresponding saturations*
    %
    %   Returned from function since since they must anyway be
    %   calculated later on

%{
Copyright 2009-2020 SINTEF Digital, Mathematics & Cybernetics.

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

    G = disc.G;

    % Mapping from all faces to internal connections
    all2int = zeros(G.faces.num,1);
    all2int(disc.internalConn) = 1:nnz(disc.internalConn);
    ix      = all2int(faces);

    % Extract cells, transmissibilities velocities, total fluxes and
    % gravity terms
    cL = disc.N(ix,1);
    cR = disc.N(ix,2);
    T  = T(ix);
    flux = flux(faces);
    g  = cellfun(@(g) g(faces), g, 'unif', false);

    % Evaluate saturations and mobilities on each side of each face
    % at each quadrature point
    nPh = numel(mob);
    [sL, sR] = deal(zeros(size(x,1),nPh));
    [rL, rR] = deal(zeros(size(x,1),nPh));
    for phNo = 1:nPh
        sL(:, phNo) = disc.evaluateDGVariable(x, cL, state, value(sdof{phNo}));
        sR(:, phNo) = disc.evaluateDGVariable(x, cR, state, value(sdof{phNo}));
    end
    if ~isempty(rdof)
        for phNo = 1:nPh
            rL(:, phNo) = disc.evaluateDGVariable(x, cL, state, value(rdof{phNo}));
            rR(:, phNo) = disc.evaluateDGVariable(x, cR, state, value(rdof{phNo}));
        end
    end
    [xL, xR] = deal(zeros(size(x,1),1));
    if ~isempty(xdof)
        xL = disc.evaluateDGVariable(x, cL, state, value(xdof));
        xR = disc.evaluateDGVariable(x, cR, state, value(xdof));
    end
    s = [sL; sR];
    sT = sum(s,2);
    s = mat2cell(s./sT, size(s,1), ones(size(s,2),1));
    r = [rL; rR];
    x = [xL; xR];
    
    for phNo = 1:nPh
        mob{phNo} = mob{phNo}([cL; cR], s, r(:,phNo), x);
    end

    % Make fake faceUpstr function
    N = [1:numel(ix); numel(ix)+1:2*numel(ix)]';
    upw = @(flag, x)faceUpstr(flag, x, N, [size(N,1), max(max(N))]);

    % Use standard MRST function to compute upstream flags
    [flagV, flagG] = getSaturationUpwind(disc.upwindType, state, g, flux, T, mob, upw);

    % For each phase, assign upstram cell and corresponding
    % saturation for each quadrature point of each face
    [upCellsV, upCellsG] = deal(repmat(cR, 1, nPh));    
    [s_v, s_G] = deal(cell(1, nPh));
    [s_v{:},s_G{:}] = deal(sR);
    for phNo = 1:nPh
        % Viscous upstream cell
        upCellsV(flagV(:,phNo), phNo) = cL(flagV(:,phNo));
        % Gravity upstram cell
        upCellsG(flagG(:,phNo), phNo) = cL(flagG(:,phNo));
    end

end
% 
% function mob = evaluateMobilities(model, mob, x, cL, cR, sdof, rsdof, rvdof)
%     
%     nPh = numel(mob);
%     for phNo = 1:nPh
%         sL(:, phNo) = disc.evaluateDGVariable(x, cL, state, double(sdof{phNo}));
%         sR(:, phNo) = disc.evaluateDGVariable(x, cR, state, double(sdof{phNo}));
%     end
%     s = [sL;sR];
%     sT = sum(s, 2);
%     if model.vapoil || model.disgas
%         rsL = disc.evaluateDGVariable(x, cL, state, double(rsdof));
%         rsR = disc.evaluateDGVariable(x, cR, state, double(rsdof));
%         rvL = disc.evaluateDGVariable(x, cL, state, double(rvdof));
%         rvR = disc.evaluateDGVariable(x, cR, state, double(rvdof));
%         mob{1} = 
%     else
%         for phNo = 1:nPh
%             mob{phNo} = mob{phNo}(s(:,phNo), sT, [cL; cR]);
%         end
%     end
%     
% end
