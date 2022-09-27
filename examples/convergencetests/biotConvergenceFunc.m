function output = biotConvergenceFunc(params, varargin)
% params is a struct with the fields:
%    Nd       : Spatial dimension (Nd = 2 or 3)
%    nref     : number of refinement level
%    eta      : Value used to set the position of the continuity point
%    gridtype : Different grid type can be used, see below
%               1 : Cartesian grid
%               2 : Triangular grid, 90 degree angles
%               3 : Equilateral triangles
%    mu       : Lame first coefficient
%    lambda   : Lame second coefficcient
%    K        : isotropic permeability, value of unique diagonal coefficient
%    alpha    : Biot coefficient
%    rho      : fluid weak compressibility coefficient

%{
Copyright 2020 University of Bergen and SINTEF Digital, Mathematics & Cybernetics.

This file is part of the MPSA-W module for the MATLAB Reservoir Simulation Toolbox (MRST).

The MPSA-W module is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

The MPSA-W module is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with the MPSA-W module.  If not, see <http://www.gnu.org/licenses/>.
%}


    opt = struct( 'blocksize', [], ...
                  'bcetazero', false);
    opt = merge_options(opt, varargin{:});
    
    Nd = params.Nd;
    nref = params.nref;
    gridtype = params.gridtype;
    eta = params.eta;
   
    output = analyticalBiot(Nd, params);
    
    % add analytical function to the params struct
    params.u_fun      = output.u_fun;
    params.p_fun      = output.p_fun;
    params.stress_fun = output.stress_fun;
    params.force_fun  = output.force_fun;
    params.src_fun    = output.src_fun;
    clear output
    
    for iter1 = 1 : nref
        
        disp(['Refinement ' num2str(iter1)])
        
        Nx = 2^iter1*ones(1, Nd); 
        G = gridForConvTest(Nx, gridtype); 

        Nc = G.cells.num;
        
        output = runBiotConvSim(G, params, 'bcetazero', opt.bcetazero, 'blocksize', opt.blocksize);
        
        % Computed displacement at cell centers
        unum = output.u;
        % Computed pressure at cell centers
        pnum = output.p;
        clear output
        
        % Analytical solution : displacement at cell centers
        uex = NaN(Nc, Nd);
        for idim = 1 : Nd
            cc{idim} = G.cells.centroids(:, idim);
        end
        
        for idim = 1 : Nd
            uex(:, idim) = params.u_fun{idim}(cc{:}); 
        end
        pex = params.p_fun(cc{:});
        
        % Errors in L2 and max - norm
        vols = G.cells.volumes;
        uerrL2(iter1) = sqrt(sum(sum(bsxfun(@times, vols.^2, (uex - unum).^2)))) / sqrt(sum(sum(bsxfun(@times, vols.^2, ...
                                                          (uex).^2))));
        perrL2(iter1) = sqrt(sum(sum(bsxfun(@times, vols.^2, (pex - pnum).^2)))) / sqrt(sum(sum(bsxfun(@times, vols.^2, ...
                                                          (pex).^2))));
        
        doplot = false;
        if doplot
            figure
            plotCellData(G, pex - pnum);
            titlename = sprintf('Difference, ref level = %d', iter1);
            title(titlename);
            colorbar
        end
        
    end
    
    % Displacement for analytical solution (last refinement)
    output.uex = uex; 
    % Pressure for analytical solution (last refinement)
    output.pex = pex;
    
    % Displacement for mpsa-mpfa solution (last refinement)
    output.unum = unum;
    % Pressure for mpsa-mpfa solution (last refinement)
    output.pnum = pnum;

    % Grid (last refinement)
    output.G = G;
    
    % L2 error between analytical and mpsa
    output.uerrL2 = uerrL2;
    output.perrL2 = perrL2;
end
