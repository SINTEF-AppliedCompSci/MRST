function output = mpsaPaperConvergenceFunc(params, varargin)
% params is a struct with the fields:
%    Nd       : spatial dimension
%    nref     : number of refinement level
%    kappa    : coefficient used for the top-corner (kappa = 1 corresponds to
%               homogeneous case)
%    alpha    : Coefficient defining lambda from mu, lambda = alpha*mu;
%    eta      : Value used to set the position of the continuity point
%    gridtype : Different grid type can be used, see below
%               1 : Cartesian grid
%               2 : Triangular grid, 90 degree angles
%               3 : Equilateral triangles
    
% Convergence test case
% title={Finite volume methods for elasticity with weak symmetry},
% author={Keilegavlen, Eirik and Nordbotten, Jan Martin},
% journal={International Journal for Numerical Methods in Engineering},
% volume={112},
% number={8},
% pages={939--962},
% year={2017},
% publisher={Wiley Online Library}

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


    opt = struct('verbose'  , false, ...
                 'blocksize', []   , ...
                 'bcetazero', false);
    
    opt = merge_options(opt, varargin{:});
    
    Nd = params.Nd;
    nref = params.nref;
    kappa = params.kappa;
    alpha = params.alpha;
    gridtype = params.gridtype;
    eta = params.eta;
    
    % at the moment, we must have mu0 = 1
    mu0 = 1;
    lamb0 = params.alpha*mu0;
   
    [u_fun, force_fun, mu_fun] = analyticalReferencePaper(Nd, kappa, mu0, lamb0);
    
    params.u_fun     = u_fun;
    params.force_fun = force_fun;
    params.mu_fun    = mu_fun;
    
    for iter1 = 1 : nref
        
        disp(['Refinement ' num2str(iter1)])
        
        Nx = 2^iter1*ones(1, Nd); 
        G = gridForConvTest(Nx, gridtype); 
        G = computeVEMGeometry(G); 
        G = computeGeometryCalc(G); 
        Nd = G.griddim;
        Nc = G.cells.num;
        
        output = runConvSim(G, params, 'bcetazero', opt.bcetazero, 'blocksize', ...
                            opt.blocksize);
        
        u = output.u;
        tbls = output.tbls;
        
        dnum = u;
        
        % Analytical solution : displacement at cell centers
        dex = NaN(Nc, Nd);
        for idim = 1 : Nd
            cc{idim} = G.cells.centroids(:, idim);
        end
        
        for idim = 1 : Nd
            dex(:, idim) = u_fun{idim}(cc{:}); 
        end
        
        % Errors in L2 and max - norm
        deL2(iter1) = sqrt(sum(sum(bsxfun(@times, G.cells.volumes.^2, (dex - dnum).^2)))) / sqrt(sum(sum(bsxfun(@times, G.cells.volumes.^2, ( dex).^2)))); 
        dem(iter1) = max(max(abs(dex - dnum))); 
        
        dostresscomputation = false;
        if dostresscomputation
            san1 = sum([s11(xf( :, 1), xf( :, 2)), s21(xf( :, 1), xf( :, 2))] ...
                       .* G.faces.normals, 2); 
            san2 = sum([s12(xf( :, 1), xf( :, 2)), s22(xf( :, 1), xf( :, 2))] ...
                       .* G.faces.normals, 2)
            stress = md.stress*reshape(dnum', [], 1); 
            
            s_ex = reshape([san1, san2]', [], 1); 
            
            sem(iter1) = max(abs(s_ex - stress))/ max(abs(stress)); 
            fa = reshape(repmat(G.faces.areas, 1, Nd)', [], 1); 
            seL2(iter1) = sqrt(sum(fa.^2 .* (stress - s_ex).^2)) / sqrt(sum(fa.^2 .* s_ex.^2)); 
        end
        
    end

    
    % displacement for analytical solution (last refinement)
    output.dex = dex; 
    % displacement for mpsa solution (last refinement)
    output.dnum = dnum;
    % L2 error between analytical and mpsa
    output.deL2 = deL2;
    
end

