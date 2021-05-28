%%

pth = mrstPath('mrst-solvent');
savepng = @(name) print([pth, '/presentation/figures/misc/', name], '-dpng', '-r300');
% savepng = @(name) [];

%%

close all

fluid = initSimpleADIFluid('n'     , [2, 2, 2], ...
                           'rho'   , [1000, 800, 100]*kilogram/meter^3, ...
                           'phases', 'WOG', ...
                           'mu'    , [1, 10, 2]*centi*poise);

sOres_i= 0.4;
fluid = addSolventProperties(fluid, 'n', 2, ...
                                    'rho', 100*kilogram/meter^3, ...
                                    'mixPar', 2/3, ...
                                    'mu'    , 1*centi*poise, ...
                                    'sOres_i', sOres_i, ...
                                     'sOres_m', 0.0);%, ...
%                                     'sSGres_i', sOres_i, ...
%                                     'sSGres_m', 0.0);


n = 100;
[sO,sS] = meshgrid(linspace(0,1,n)', linspace(0,1,n)');

so = sO;
ss = sS;
% too_much = sO + sS > 1;
% sO(too_much) = 0;
% sS(too_much) = 0;

% sS = 1 - sO;
% sG = 0.0*ones(size(sS));
sW=0.0*ones(size(sS));
nd=size(sS);
p=0*barsa*ones(size(sS));
p = reshape(p, [], 1);
%sG+sW+sO+sS=1
sG = 1-(sW+sO+sS);
sg = sG;
%%{
str={'W','O','G','S'};
for i=1:3
    sat{i}=reshape(eval(['s',str{i}]),[],1);
end
[sW,sO,sG] = initVariablesADI(sat{1},sat{2},sat{3});
sS=1-(sW+sO+sG);
%}

mobMult = 1;
[sWres, sOres , sSGres ] = computeResidualSaturations(fluid, p , sG , sS );
[krW, krO, krG, krS] = computeRelPermSolvent(fluid, p, sW, sO, sG, sS, sWres, sOres, sSGres, mobMult);
[muW, muO, muG, muS, rhoW, rhoO , rhoG , rhoS ] ...
    = computeViscositiesAndDensities(fluid, p , sO , sG , sS , sOres , sSGres );
[bW , bO , bG , bS ] = computeFormationVolumeFactors(fluid, p , rhoO , rhoG , rhoS );
str={'O','G','S'};
mobcell = cell(3,4);
% sy = reshape(double(sG), nd);


ssp = deal(linspace(0,1,100)');
b = sOres_i + 2*ssp - 1;
sgp = (-b + sqrt(b.^2 - 4*(ssp.^2-ssp)))/2;
sor = sgp./(ssp + sgp).*sOres_i;

[krw, kro, krg, krs] = computeRelPermSolvent(fluid, 0, 0, sor, sgp, ssp, 0, sor, 0, 1);
kro = kro + 1e-2;
krg = krg + 0;
krs = krs + 0;

% t = delaunay(sx, sy);



azel = repmat([-40, 30], 3, 1);
azel(2,:) = azel(2,:) + [180,0];

for i=1:3
    kr=eval(['kr',str{i}]);
    mu=eval(['mu',str{i}]);
    b=eval(['b',str{i}]);
    mob=kr./mu;
    mobb=b.*mob;
    
%     mobbv=reshape(double(mobb),nd);    
    mobbv=reshape(double(kr),nd);
%     mobbv=reshape(double(mu),nd);
%     mobbv=reshape(double(b),nd);
   

    if 0
        figure(i)
        subplot(1,4,1),cla    
    else
        figure
        [mapx, mapy] = ternaryAxis('names', {'S_g', 'S_s', 'S_o'});
    end
    
    mobbv(sg < 0) = nan;
%     trisurf(t, sx, sy, mobbv);
    contourf(mapx(sg, ss, so), mapy(sg,ss,so), mobbv, 20, 'linecolor', 0.5.*[1,1,1])
    hold on
    [mapx, mapy] = ternaryAxis('names', {'S_g', 'S_s', 'S_o'});
    plot(mapx(sgp, ssp, sor), mapy(sgp,ssp,sor), 'color', 0.99*[1,1,1], 'linewidth', 2)
    xlabel('S_o');
    ylabel('S_s');
%     zlabel(['k_{r', lower(str{i}), '}']);
    hold on
    kro = eval(['kr',lower(str{i})]);
%     plot3(sor, ssp, kro, 'r', 'linewidth', 2);
    hold off
%     axis([0 1 0 1]);
    axis equal tight
    caxis([0,1])
    colormap(parula)
    if 1
        ax = gca;
%         ax.ZTick = [0,0.5,1];
        ax.XLabel.FontSize = 20;
        ax.YLabel.FontSize = 20;
%         ax.ZLabel.FontSize = 20;
        savepng(['kr', str{i}]);
    end
        
    
    mobcell{i,1} = mobbv;
    
    nn = zeros(4,1);
    
    nn(1) = nnz(isnan(mobbv));
    for k=1:3
        
        
        dd=diag(mobb.jac{k});
%         dd = diag(kr.jac{k});
%         dd=diag(mu.jac{k});
%         dd=diag(b.jac{k});
        
        nn(k+1) = nnz(isnan(full(dd)));
        mobbd=reshape(dd,nd);
        if 0
            subplot(1,4,1+k),cla
            surf(sg, ss, mobbd)
        end
        
        mobcell{i,k+1} = mobbd;
        
    end    
    
    
    fprintf('%s:', str{i})
    for i = 1:4
        fprintf('\t%d', nn(i));
    end
    fprintf('\n');
    
end

%%
% <html>
% <p><font size="-1">
% Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.
% </font></p>
% <p><font size="-1">
% This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).
% </font></p>
% <p><font size="-1">
% MRST is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% </font></p>
% <p><font size="-1">
% MRST is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% </font></p>
% <p><font size="-1">
% You should have received a copy of the GNU General Public License
% along with MRST.  If not, see
% <a href="http://www.gnu.org/licenses/">http://www.gnu.org/licenses</a>.
% </font></p>
% </html>
