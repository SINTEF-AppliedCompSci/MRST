mrstModule add vemmech mpsaw vem mpfa

close all

params = struct('nref'    , 1, ...
                'Nd'      , 2, ...
                'gridtype', 1, ...
                'eta'     , 0, ...
                'mu'      , 1, ...
                'lambda'  , 1, ...
                'alpha'   , 1, ...
                'K'       , 1, ...
                'tau'     , 1, ...
                'rho'     , 0);

% Compute numerical solution
output1 = biotConvergenceFunc(params);
output2 = biotConvergenceFunc(params, 'blocksize', 2);

fprintf('Displacement error (standard) : %g\n', output1.uerrL2);
fprintf('Displacement error (block assembly): %g\n', output2.uerrL2);

return

u = output.u;
p = output.p;

dotest = true;
if dotest
    tbls = output.tbls;
    isBoundary = any(G.faces.neighbors == 0, 2); 
    bcfaces =  find(isBoundary);
    bcfacetbl.faces = bcfaces;
    bcfacetbl = IndexArray(bcfacetbl);
    cellfacetbl = tbls.cellfacetbl;
    bccellfacetbl = crossIndexArray(cellfacetbl, bcfacetbl, {'faces'});
    bccelltbl = projIndexArray(bccellfacetbl, {'cells'});
    bccells = bccelltbl.get('cells');
    bcp = p(bccells);
end    

doplot = true;
if doplot
    
    figure
    plotCellData(G, p);
    title('Pressure, numerical solution');
    
    figure
    plotCellData(G, u(:, 1));
    title('Displacement, x-direction, numerical solution');
    
    figure
    plotCellData(G, u(:, 2));
    title('Displacement, y-direction, numerical solution');

    % prepare input for analytical functions
    for idim = 1 : d
        cc{idim} = G.cells.centroids(:, idim);
    end
    
    figure
    p_exact = params.p_fun(cc{:});
    plotCellData(G, p_exact);
    title('Pressure, analytical solution');
    
    figure
    u1_exact = params.u_fun{1}(cc{:});
    plotCellData(G, u1_exact);
    title('Displacement, x-directiom, analytical solution');

    figure
    u2_exact = params.u_fun{2}(cc{:});
    plotCellData(G, u2_exact);
    title('Displacement, y-directiom, analytical solution');

end

%% Copyright Notice
%
% <html>
% <p><font size="-1">
% Copyright 2020 University of Bergen and SINTEF Digital, Mathematics & Cybernetics.
% </font></p>
% <p><font size="-1">
% This file is part of the MPSA-W module for the MATLAB Reservoir Simulation Toolbox (MRST).
% </font></p>
% <p><font size="-1">
% The MPSA-W module is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% </font></p>
% <p><font size="-1">
% The MPSA-W module is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% </font></p>
% <p><font size="-1">
% You should have received a copy of the GNU General Public License
% along with the MPSA-W module.  If not, see
% <a href="http://www.gnu.org/licenses/">http://www.gnu.org/licenses</a>.
% </font></p>
% </html>

