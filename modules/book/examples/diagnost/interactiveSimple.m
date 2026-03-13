%% Simple 2D reservoir with two injectors and three producers
% This is an example discussed in the MRST book, which is a continuation
% of the 'showDiagnostBasics' example, in which we have added two
% low-porosity regions as a simple representation of two sealing faults,
% moved the two injectors slightly, and switched all wells from rate to
% bottom-hole pressure control.
%
% In the following, we give instructions for how you can use the diagnostic
% tool to manipulate the well controls to improve the overall sweep of the
% reservoir.

%% Set up the geomodel and specify wells
[nx,ny] = deal(64);
G = cartGrid([nx,ny,1],[5000,2500,10]);
G = computeGeometry(G);
p = gaussianField(G.cartDims(1:2), [0.2 0.4], [11 3], 2.5);
p(round(2*end/3):end,round(end/3)) = 1e-3;
p(1:round(end/3),round(2*end/3)) = 1e-3;
K = p.^3.*(1.5e-5)^2./(0.81*72*(1-p).^2);
rock = makeRock(G, K(:), p(:));

%% Set up and solve flow problem, compute diagnostics
n = 12;
W = addWell([],  G, rock, nx*(n-6)+n/2+1, ...
    'Type', 'bhp', 'Comp_i', 1, 'name', 'I1', 'Val', 200*barsa);
W = addWell(W, G, rock, nx*(n-6)+n/2+1+nx-n, ...
    'Type','bhp',  'Comp_i', 1, 'name', 'I2', 'Val', 200*barsa);
W = addWell(W, G, rock, round(G.cells.num-(n/2-.5)*nx-.3*nx), ...
    'Type','bhp',  'Comp_i', 0, 'name', 'P1', 'Val', 100*barsa);
W = addWell(W, G, rock, G.cells.num-(n-.5)*nx, ...
    'Type','bhp',  'Comp_i', 0, 'name', 'P2', 'Val', 100*barsa);
W = addWell(W, G, rock, round(G.cells.num-(n/2-.5)*nx+.3*nx), ...
    'Type','bhp',  'Comp_i', 0, 'name', 'P3', 'Val', 100*barsa);

%%
close all;
mrstModule add diagnostics
interactiveDiagnostics(G, rock, W, 'showGrid', true);
axis normal tight; view(0,90);

%% Plot of reservoir model:
%  - select the displayed property to be 'porosity'
%  - select the displayed property to be 'drainage blend'
% To add a fancy colorbar to the porosity plot, you can use the following
% call:
%  [hc,hh]=colorbarHist(rock.poro,[min(rock.poro) max(rock.poro)],'South');
% and then manipulate the 'position' property of each handle to improve the
% placement of the colorbar. To confirm the flux allocation for each of the
% wells, you simply click on the text that signifies each of the wells.

%% Plot of time evolution
% - Set display to be 'forward TOF', and selection to be 'Flood volumes'
% - set 'Max TOF' to be 25, 75, 125, and 175
% - or use the 'Play TOF' button to play the evolution of neutral
%   displacement fronts

%% Manual optimization of flow pattern
% First, change the colormap to
colormap(gray.^10);
% and then switch the displayed property to be 'sum of tofs'

% From the 'plots' tab, you can first plot the F-Phi diagram, and then say
figure(3); hold on
% to ensure that if you push this button subsequently, the new F-Phi
% diagram will be added on top of the preexisting one(s) rather than
% replacing it. (Depending upon what you have done before you come to this
% point, the number of the figure may be different. Check the window header
% to get the right number.)

% Case 1:
% Use the 'edit wells' button from the 'plots' tab to edit the various
% wells. We can start by increasing the pressure of 'P2' from 100 bar to
% 130 bar to decrease the flow rate in the I1-P2 region. Once you have
% entered the new value (13000000), you push 'apply' and a new flow field
% with flow diagnostics will be computed. You can now change the displayed
% property to be 'sum of TOFs', and use the 'F/Phi diagram' button from the
% 'plots' tab to plot the F-Phi diagram and compute the Lorenz coefficient
% of the new well setup. (To distinguish the curves, you may want to
% manually change the color of the new graph using the 'Edit plot' tool
% from the plotting window.)

% Case 2:
% To increase the flow in the I2-P3 region, we can repeat the exercise and
% set P3 to operate at 80 bar.

% Case 3:
% To also sweep the large unswept region east of P3, we can introduce a new
% well in the south-east of this region, just north of the sealing fault.
% To do this, we use the 'add' button from the 'edit wells' workflow. This
% will bring up a new window in which you can click the cell where you want
% to place the well. When you click a cell, a new window will pop up; this 
% can be ignored. Once you are happy with the placement, you close the
% editing window, and then you can sett the correct value in the new well.
% In the book, we propose to set the new well to operate at 140 bar, and
% similarly increase the pressure of I2 to 220 bar. Then you push 'apply',
% go back to view the 'sum of TOFs' and produce a new F-Phi diagram

%% Plot of time evolution
% - Set display to be 'forward TOF', and selection to be 'Flood volumes'
% - set 'Max TOF' to be 25, 75, 125, and 175
% - or use the 'Play TOF' button to play the evolution of neutral
%   displacement fronts

%% Exercise:
% Can you make any suggestions for further improvements?
