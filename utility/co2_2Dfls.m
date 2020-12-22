% Function for the co2 assesment of the 2D leakage problem.
%
%{
Copyright 2020, NORCE Norwegian Research Centre AS, Computational 
Geosciences and Modeling. 

This file is part of the ad-micp module.

ad-micp is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

ad-micp is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this file.  If not, see <http://www.gnu.org/licenses/>.
%}
function [state,lr,tt,Q1]=co2_2Dfls(statesMICP) 

    % Setup some schedule
    dt = .5;
    nt = 120+59+23+24*(100-1); % 100 days
    clear schedule
    timesteps = repmat(dt, nt, 1);

    % Grid
    L = 500;
    H = 160;
    [X1,Y1] = meshgrid([-L:10:-50 180:10:L], [0:5:30 31:1:129 130:5:H]);
    [X2,Y2] = meshgrid(-50:10:50, 0:5:30);
    [x,y] = meshgrid([100-.3 100 100+.3], 0:.25:H);
    [xc,yc] = meshgrid([-1 0 1], 130:.25:H);
    [xw,yw] = meshgrid([-50*exp(0:-0.25:-3.6) 50*exp(-3.6:0.25:0)], ...
                                                                  130:1:H);
    [xl,yl] = meshgrid([100-50*exp(0:-0.25:-5) 100+80*exp(-5:0.125:0)], ...
                                                                    0:1:H);
    [xwc1,ywc1] = meshgrid(50*exp(-3.6:0.25:0), 130:.25:137.5);
    [xwc2,ywc2] = meshgrid(50*exp(-3.6:0.25:0), 137.5:.5:145);
    [xwc3,ywc3] = meshgrid(100-50*exp(0:-0.25:-5), 130:.25:137.5);
    [xwc4,ywc4] = meshgrid(100-50*exp(0:-0.25:-5), 137.5:.5:145);
    P = unique([X1(:) Y1(:); X2(:) Y2(:); x(:) y(:); xw(:) yw(:); xl(:) ...
          yl(:); xc(:) yc(:); xwc1(:) ywc1(:); xwc2(:) ywc2(:); xwc3(:) ...
                                        ywc3(:); xwc4(:) ywc4(:)], 'rows');
    G = triangleGrid(P);
    G = computeGeometry(G);
    c = G.cells.centroids;
    G = removeCells(G, (c(:,1)<100-.3 | c(:,1)>100+.3) & (c(:,2)<130 & ...
                                                               c(:,2)>30));
    G = makeLayeredGrid(pebi(G),1);
    G = computeGeometry(G);
    c = G.cells.centroids;
    G = removeCells(G, (c(:,1)<99.9 | c(:,1)>100.1) & (c(:,2)<130 & ...
                                                               c(:,2)>30));
    G = computeGeometry(G);
    c = G.cells.centroids;
    G = removeCells(G, c(:,1)<-.5-eps & c(:,2)>130);
    G = computeGeometry(G);
    c = G.cells.centroids;
    G = removeCells(G, c(:,1)<99.8 & c(:,2)<30);
    G = computeGeometry(G);
    C = ones(G.cells.num,1);

    % Rock
    K0 = 2e-14;
    K = K0.*ones(G.cells.num,1);
    cellsfrac =  G.cells.indexMap;
    cellsfrac1 = cellsfrac(G.cells.centroids(:,1)>99.9 & ...
            G.cells.centroids(:,1)<100.1 & G.cells.centroids(:,2)<130 & ...
                                                G.cells.centroids(:,2)>30);
    cellsF =  G.cells.indexMap;
    idx = ismember(cellsF,cellsfrac1);
    K(idx) = 1e-12;
    K0=K;
    porosity = 0.15;

    % Model properties
    fluid.cells = C;
    fluid.eta = 3;
    fluid.crit = .1;
    fluid.kmin = 1e-20;
    fluid.K = @(poro) (K0.*((poro-fluid.crit)/(porosity-fluid.crit))...
        .^fluid.eta+fluid.kmin).*K0./(K0+fluid.kmin).*(poro>fluid.crit)+...
                                            fluid.kmin.*(poro<=fluid.crit);
    
    % Compute current porosity and permeability                                    
    poro = 0.15-statesMICP.c-statesMICP.b;
    KK = fluid.K(poro);
    rock = makeRock(G, KK, poro);

    % Fluid parameters
    fluid.muW   = 2.535e-4;       % Water viscosity 
    fluid.muO   = 3.95e-5;        % CO2 viscosity                                       
    fluid.bW    =  @(p) 0*p + 1;  % Water formation volume factor
    fluid.bO    =  @(p) 0*p + 1;  % CO2 formation volume factor 
    fluid.rhoWS = 1045;           % Water density
    fluid.rhoOS = 479;            % CO2 density  

    % Gravity
    gravity on
    gravity y

    % Create Well
    Q1 = (1600/day)/L; %m^3/day
    cellsWell =  1:1:G.cells.num;
    cellsWell1 = cellsWell(G.cells.centroids(:,1)>-.5 & ...
                   G.cells.centroids(:,1)<.5 & G.cells.centroids(:,2)>130);
    W = addWellMICP([], G, rock, cellsWell1, 'Type', 'rate', 'Comp_i', ...
         [eps,1-eps], 'Val', Q1, 'Radius',.15*meter,'name', 'I','dir','y');

    % Create model
    model = CO2Model(G, rock, fluid);

    % Boundary condition
    f = boundaryFaces(G);
    f = f(abs(G.faces.normals(f,1))>eps & (G.faces.centroids(f,1)<-L+2 ...
                                           | G.faces.centroids(f,1)>L-2 ));
    fp = G.faces.centroids(f,2)*fluid.rhoWS*norm(gravity);
    bc = addBC([], f, 'pressure', fp, 'sat', [0 0]);
    
    % Well different rates and times
    N = 2; % Number of injection changes
    M = zeros(N,2); % Matrix where entries per row are: time and rate.
    M(1,1) = 121*dt; 
    M(1,2) = Q1;
    M(2,1) = (121+59)*dt; 
    M(2,2) = Q1;
    
    % Make schedule
    schedule = simpleSchedule(timesteps,'W',W,'bc',bc);
    for i=1:N
        schedule.control(i+1)=schedule.control(i);
        schedule.control(i+1).W(1).val=M(i,2);
        schedule.control(i+1).W(1).sign=sign(M(i,2));
        schedule.step.control(M(i,1)/dt:end)=i+1;
        schedule.step.val(M(i,1)/dt:end)=minute;
    end
    schedule.step.val(M(2,1)/dt:end)=hour;

    % Initial state
    state0 = initState(G, W, G.cells.centroids(:,2)*fluid.rhoWS* ...
                                      norm(gravity), [1-eps, eps]);

    % Simulate
    [~, states] = simulateScheduleAD(state0, model, schedule);
    state=states{end};

    % Compute leakage rate
    cellsfa =  1:1:G.faces.num;
    cellsfac = cellsfa(G.faces.centroids(:,2)<80.45+1000*eps & ...
                                    G.faces.centroids(:,2)>80.35-1000*eps);
    for i=1:1:nt
        lr(i)=abs(states{i}.flux(cellsfac(1),2));
    end
    for i=1:1:120
        tt(1,i)=i*.5/day;
    end
    for i=121:1:120+59
        tt(1,i)=120*.5/day+(i-120)*minute/day;
    end
    for i=120+60:1:nt
        tt(1,i)=120*.5/day+59*minute/day+(i-120-59)*hour/day;
    end
end