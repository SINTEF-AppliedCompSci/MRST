function plotSolventFluidProps(model, propNames, phases, varargin)
% Terneary plots of HC components in black-oil solvent model.

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

    opt = struct('np', 100, 'pressure', 100*barsa);
    
    opt  = merge_options(opt, varargin{:});
    
    p = opt.pressure;

    np = opt.np;
    
    df = get(0, 'defaultFigurePosition');
    fluid = model.fluid;
    
    [sO,sS] = meshgrid(linspace(0,1,np)', linspace(0,1,np)');
    sG = 1-(sO+sS); sW = sO.*0;
    
    [sor, sgc] = solveSr(fluid,p);
    ss  = linspace(0,1,100)';
    so = 1 - ss - sgc;
    sg = 1 - ss - sor;
    
    [rs, rv] = deal(zeros(np*np,1));
    sat = true(np*np,1);
    
    sW = sW(:); sO = sO(:); sG = sG(:); sS = sS(:);
    [sWr, sOr, sGc]      = computeResidualSaturations(model, p, sW, sG, sS);
    [krW, krO, krG, krS] = computeRelPermSolvent(model, p, sW, sO, sG, sS, sWr, sOr, sGc, 1);
    [muW, muO, muG, muS, rhoW, rhoO, rhoG, rhoS, bW, bO, bG, bS, pW, pG] = computeViscositiesAndDensities(model, p.*ones(np*np,1), sW, sO , sG , sS , sOr , sGc, rs, rv, sat, sat);
    sO = reshape(sO, np,np); sG = reshape(sG, np,np); sS = reshape(sS, np,np);

    for propNo = 1:numel(propNames)
        
        figure('Position', [df(1:2), 1000, 400]);

        for phNo = 1:numel(phases)

            subplot(1,3,phNo)

            name = [propNames{propNo}, phases{phNo}];
            prop = reshape(eval(name),[np,np]);
            prop(sG<0) = nan;

            [mapx, mapy] = ternaryAxis('names', {'S_g', 'S_s', 'S_o'});
            contourf(mapx(sG, sS, sO), mapy(sG,sS,sO), prop, 20, 'linecolor', 0.5.*[1,1,1])
            [mapx, mapy] = ternaryAxis('names', {'S_g', 'S_s', 'S_o'});
            plot(mapx(sg(sg>=0), ss(sg>=0), sor(sg>=0)), mapy(sg(sg>=0),ss(sg>=0),sor(sg>=0)), 'color', 0.99*[1,1,1], 'linewidth', 2)
            plot(mapx(sgc(so>=0), ss(so>=0), so(so>=0)), mapy(sgc(so>=0),ss(so>=0),so(so>=0)), 'color', 0.99*[1,1,1], 'linewidth', 2)
            
            axis([0,1,0,sqrt(1-0.5^2)]); axis equal
            colorbar('location', 'southOutside');
            title(name, 'position', [0.5,-0.5]);   

        end

    end

end

function [sor, sgc] = solveSr(fluid, p)

    n = 100;

    ss = linspace(0, 1, n)';
    
    M = @(sor) fluid.Ms(fluid.satFrac(ss,1-sor)).*fluid.Mp(p);

    fun = @(sor) fluid.sOr_i.*(1 - M(sor)) ...
                + fluid.sOr_m(0).*M(sor) - sor;
    
    sor = initVariablesADI(fluid.sOr_i*ones(n,1));
    
    tol = 1e-10;
    eqs = fun(sor);
    err = norm(value(eqs), inf);
    while err > tol
        
        dsr = -eqs.jac{1}\eqs.val;
        sor = sor + dsr;
        
        eqs = fun(sor);
        err = norm(value(eqs), inf);
        
    end
    
    M = @(sgc) fluid.Ms(fluid.satFrac(ss,ss + sgc)).*fluid.Mp(p);

    fun = @(sgc) fluid.sGc_i.*(1 - M(sgc)) ...
                + fluid.sGc_m(0).*M(sgc) - sgc;
    
    sgc = initVariablesADI(fluid.sGc_i*ones(n,1));
    
    eqs = fun(sgc);
    err = norm(value(eqs), inf);
    while err > tol
        
        dsc = -eqs.jac{1}\eqs.val;
        sgc = sgc + dsc;
        
        eqs = fun(sgc);
        err = norm(value(eqs), inf);
        
    end
    sor = value(sor);
    sgc = value(sgc);
end
                 
                 
                 
                 