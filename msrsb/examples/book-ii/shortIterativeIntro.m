%% Short introduction to iterative multiscale method
% This example continues from the shortMultiscaleIntro.m script and
% contains the code necessary to produce the plots presented in Section 2.6
% of the multiscale tutorial chapter.
if exist('ranShortMultiscaleIntro','var')~=1
    shortMultiscaleIntro;
end

%% Set up global multiscale stage and smoother stage
[A,b]    = deal(state.A, state.rhs);
mssolve  = @(d) basis.B * mldivide(basis.R*A*basis.B, basis.R*d);
fn       = getSmootherFunction('type', 'ilu0');
S        = fn(A, b);
resnorm  = @(p) norm(b-A*p,2)/norm(b,2);

%% Plot of residuals
figure('Position',[300 340 1000 420]);
titles = {'First multiscale','After smoothing','Second multicale','After smoothing'};
for i=1:4
    subplot(2,2,i); cla
    switch i
        case 1    
            p = mssolve(b);
        case {2,4} 
            p = p + S(b-A*p);
        case 3
            p = p + mssolve(b-A*p);
    end
    plotCellData(G, abs(A*p-b),'EdgeColor','none'); axis tight
    plotFaces(CG,1:CG.faces.num,'linewidth',1);
    axis equal tight, set(gca,'FontSize',12);
    caxis([0 2e-8]);
    title(titles{i},'FontWeight','normal');
end
colormap(.65*flipud(parula.^6)+.35)
set(colorbar,'Position',[.9 .11 .025 .33]);

%% Two-stage version compared with GMRES
nit=65;
res = zeros(nit+1,1);
p = mssolve(b);
res(1) = resnorm(p);
for i=2:nit+1
    p = p + S(b - A*p);
    p = p + mssolve(b-A*p);
    res(i) = resnorm(p);
end

% Use GMRES
[pgmres, rep] = solveMultiscaleIteratively(A, b, [], basis, fn, 5e-15, nit, @mldivide, true);
%[state_it,rep] = incompMultiscale(state0, CG, hT, fluid, basis, ...
%    'bc', bc,'getSmoother',fn,'iterations',nit, 'useGMRES', true, 'tolerance',5e-15);

figure
plot(1:nit+1,res,'-o',1:nit+1,rep.resvec,'-s','MarkerFaceColor',[.8 .8 .8]); 
set(gca,'YScale','log'); legend('MsRSB+ILU(0)', 'MsRSB+ILU(0), GMRES');
xlabel('Iteration number'),ylabel('Normalized residual');

%%
% <html>
% <p><font size="-1">
% Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.
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
