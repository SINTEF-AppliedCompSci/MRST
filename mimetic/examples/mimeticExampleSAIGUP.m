%% Effects of Transmissibility Multipliers
% Consider once more the incompressible, immiscible two-phase pressure
% equation
%
% $$\nabla\cdot v = q, \qquad v=\textbf{--}\lambda K\nabla p$$
%
% One approach to modelling conduits and barriers to flow is to introduce
% scalar multipliers on the transmissiblity of a set of grid block
% connections (cell interfaces).  While historically strongly tied to the
% application of a particular discretisation method (the two-point flux
% approximation), recent development has demonstrated how to incorporate
% such multipliers into other consistent and convergent discretisation
% schemes such as the mimtic method.
%
% Here we illustrate effects of transmissibility mulipliers on the
% qualitative behaviour of the solution of the above problem.  Moreover, we
% demonstrate how the practitioner may introduce multiplier effects in the
% mimetic method as implemented in MRST.
%
% <html>
% We use a realisation from the <a
% href="http://www.nr.no/pages/sand/area_res_char_saigup">SAIGUP</a> study
% and continue the <a href="../../grids/html/saigupModelExample.html">grid
% example</a> from a previous demonstration.
% </html>

mrstModule add mimetic incomp deckformat
mrstVerbose false
gravity off

%% Read input file and construct model
% The SAIGUP model is one of the standard data sets provided with MRST. We
% therefore check if it is present, and if not, download and install it.
% Then, we read the ECLIPSE input and construct the necessary data objects to
% make a simulation model

fn = fullfile(getDatasetPath('SAIGUP'),'SAIGUP.GRDECL');
grdecl = readGRDECL(fn);
grdecl = convertInputUnits(grdecl, getUnitSystem('METRIC'));
G      = processGRDECL(grdecl);
G      = computeGeometry(G);
rock   = grdecl2Rock(grdecl, G.cells.indexMap);

%% Modify the permeability to avoid singular tensors
% MRST, or rather the mimetic method implemented in MRST, requires a
% positive definite permeability tensor in all active grid blocks.  The
% input data of the SAIGUP realisation, however, has zero vertical
% permeability in a number of cells.  We work around this issue by
% (arbitrarily) assigning the minimum positive vertical (cross-layer)
% permeability to the grid blocks that have zero cross-layer permeability.
%
% We do emphasise that the above restriction does mean that the mimetic
% method currently implemented in MRST is not capable of handling fully
% sealing barriers, at least if those barriers are represented as zero
% grid block permeability.
is_pos = rock.perm(:, 3) > 0;
rock.perm(~is_pos, 3) = min(rock.perm(is_pos, 3));

%% Extract multipliers from the input data
% The <matlab:doc('computeTranMult') computeTranMult> function, new in MRST
% release 2011a, extracts multiplier data from input, typically represented
% by the ECLIPSE keywords |MULTX|, |MULTY| or |MULTZ|, and computes a
% scalar multiplier value (default |1.0|) for each reservoir connection.
% The result, |m|, is represented as a scalar for each face for each
% grid block, with internal faces being represented twice.
%
% In a sense, the |computeTranMult| function plays a rôle similar to the
% |grdecl2Rock| function that extracts and reshapes permeability data from
% the input model information.
m = computeTranMult(G, grdecl); clear grdecl

%% Inspect multiplier structure
% The multiplier data |m| is provided once for each face in each grid cell,
% meaning internal faces are represented twice.  However, we need one value
% per unique face (connection) when we're plotting this data.  Fortunately,
% the |computeTranMult| function computes a symmetric assignment, so we
% need only extract the last multiplier value pertaining to any particular
% face.  This is accomplished by indexing into the (smaller) result array
% |mface| with the cell-to-face mapping and then doing regular assignment.
mface = zeros([G.faces.num, 1]);
mface(G.cells.faces(:,1)) = m;

%%
% We display the grid structure as a background and add some transparency
% to the graphics so as not to obscure the multiplier data.  Moreover, we
% only display the multipliers that introduce barriers, i.e., the
% multipliers that are less than unity.  To accentuate the structure of
% this data, we plot the logarithm of the actual multiplier values.
newplot
[az, el] = deal(-51, 30); daspctrat = [1.2, 2, 1/2];

plotGrid (G, 'EdgeAlpha', 0.05, 'FaceColor', 'none');
plotFaces(G, find(mface < 1), log10(mface(mface < 1)));

view(az, el)
axis tight off, colorbar NorthOutside

set(gca, 'DataAspectRatio', daspctrat);
zoom(1.25)

%%
% We notice that these multipliers almost exclusively introduce partially
% sealing flow barriers in the vertical direction.  We will additionally
% introduce lateral flow barriers.  The |processGRDECL| function marks
% non-neighbouring connections associated with a non-matching corner-point
% description using a non-zero tag.  We first inspect the geometrical
% structure spanned by these non-neighbouring connections.
newplot
plotGrid (G, 'EdgeAlpha', 0.05, 'FaceColor', 'none');
plotFaces(G, find(abs(G.faces.tag) > 0), ...
          'FaceColor', 'Red', 'EdgeColor', 'none');
view(-126, 40), set(gca, 'DataAspectRatio', daspctrat);
axis tight off
zoom(1.5)

%%
% We observe that the fault structures span much of the lateral portion of
% the model.  We will introduce significant flow restriction by assigning a
% transmissibility multiplier of $5\cdot 10^{\textbf{--} 5}$ to faces
% associated with any of the faults.
fault_mult = ones([G.faces.num, 1]);
fault_mult(abs(G.faces.tag) > 0) = 5.0e-5;

fault_mult = fault_mult(G.cells.faces(:,1));

%% Define driving forces
% We define a water-flooding scenario in which we inject one 50,000th of
% the total model pore volume per day in a single injection well.  We
% similarly define a production well controlled by a bottom-hole pressure
% target of 150 bars.
tot_pv     = sum(poreVolume(G, rock));
inj_rate   = (tot_pv / 50e3) / day;
prod_press = 150*barsa;

%%
% We will discretise the pressure equation using both the mimetic and the
% TPFA method.  These methods require different empirical constants in the
% Peaceman well model, here representeded by the |'InnerProduct'| option to
% the <matlab:doc('verticalWell') verticalWell> well constructor function.
%
% First we construct the injection and production wells using the
% |'ip_quasirt'| inner product.  This inner product is applicable to the
% mimetic method we construct below.
%
% We place the injector in an area of the model that is mostly boxed in by
% the fault structure, so we expect the faults to significantly affect the
% resulting reservoir flows.  Moreover, the producer is in an area that is
% also mostly shielded from the remainder of the reservoir.
Wm = verticalWell([], G, rock, 30, 85, 1:10       , ...
                 'Type', 'rate', 'Val', inj_rate, ...
                 'Radius', 12.5*centi*meter     , ...
                 'Name', 'Inj'                  , ...
                 'InnerProduct', 'ip_quasirt', 'Comp_i', [1, 0]);

Wm = verticalWell(Wm, G, rock,  5,  12, 1:10       , ...
                 'Type', 'bhp', 'Val', prod_press, ...
                 'Radius', 12.5*centi*meter      , ...
                 'Name', 'Prod'                  , ...
                 'InnerProduct', 'ip_quasirt', 'Comp_i', [0, 1]);

plotWell(G, Wm, 'Height', 200, 'Color', 'magenta', ...
         'cylpts', 30, 'FontSize', 20);

%%
% Secondly, we construct the injection and production wells using the
% |'ip_tpf'| inner product.  This is Peaceman's traditional well model that
% contains the empirical factor |0.14| and is applicable to the two-point
% discretisation defined below.
Wtp = verticalWell([], G, rock, 30, 85, 1:10   , ...
                   'Type', 'rate', 'Val', inj_rate, ...
                   'Radius', 12.5*centi*meter     , ...
                   'Name', 'Inj'                  , ...
                   'InnerProduct', 'ip_tpf', 'Comp_i', [1, 0]);

Wtp = verticalWell(Wtp, G, rock,  5,  12, 1:10    , ...
                    'Type', 'bhp', 'Val', prod_press, ...
                    'Radius', 12.5*centi*meter      , ...
                    'Name', 'Prod'                  , ...
                    'InnerProduct', 'ip_tpf', 'Comp_i', [0, 1]);

%% Calculate face transmissibilites
% In this problem we consider only barriers to flow.  These manifest as
% transmissibility multipliers less than unity.  Selecting multipliers less
% than $10^{\textbf{--} 3}$ we pick out those multipliers whose effect on
% the overall solution will be the most distinctive.
%
% An upcoming paper in the SPE Journal shows that multipliers can be
% incorporated into the mimetic method by computing an associated face
% transmissibility and subsequently modifying the mimetic inner product
% using this transmissibility.  In particular, if the connection (face)
% between cells |i| and |j| is affected by a multiplier $m \in (0,1)$, then
% we may add a diagonal contribution of
%
% $$\mathsf{T}_{ij} \, \frac{m}{1\textbf{--} m} $$
%
% to the corresponding term in the mimetic inner product in cells |i| and
% |j|, respectively.
%
% Here we incorporate the synthetic fault transmissibilities defined above
% into the actual multiplier data defined in the input, compute the face
% transmissibilities $\mathsf{T}_{ij}$ (result |ftrans|) and the multiplier
% update outlined above (result |facetrans|). Function |computeMimeticIP|
% will then incorporate these values into the pertinent inner products.
m = m .* fault_mult;
i = m < 1e-3;
T = computeTrans(G, rock);

ftrans    = 1 ./ accumarray(G.cells.faces(:,1), 1 ./ T);
facetrans = ftrans(G.cells.faces(i,1)) .* (m(i) ./ (1 - m(i)));


%%
% We now calculate the mimetic inner product with and without explicit face
% transmissibility.  Note that we employ the |'ip_quasirt'| inner product
% when constructing the discretisations.  This is consistent with the well
% object defined above.

S0 = computeMimeticIP(G, rock, 'InnerProduct', 'ip_quasirt');
S1 = computeMimeticIP(G, rock, 'InnerProduct', 'ip_quasirt', ...
                      'FaceTrans', ...
                      [G.cells.faces(i,1), facetrans]);

%% Create a reservoir fluid
% This model problem uses a standard, incompressible fluid with a factor
% |10| mobility ratio and quadratic relative permeability curves without
% residual effects.
fluid = initSimpleFluid('mu' , [   1,  10]*centi*poise     , ...
                        'rho', [1000, 700]*kilogram/meter^3, ...
                        'n'  , [   2,   2]);

%% Initialize the reservoir state object.
% The function <matlab:doc('initState') initState> function, new in MRST
% 2011a, initialises a reservoir- and well solution state object whilst
% checking its input parameters for basic consistency.  In particular, the
% number of phases must be the same in the well objects and the reservoir
% state object.
%
% We define the reservoir initial condition as a constant pressure of 100
% bars and filled with oil (zero water saturation).
xm  = initState(G, Wm   , 100*barsa, [0, 1]);
xtp = initState(G, Wtp, 100*barsa, [0, 1]);

%% Solve flow problems with and without multiplier effects
% Mimetic discretisation *without* multipliers.
sols{1} = incompMimetic(xm, G, S0, fluid, 'Wells', Wm);

% Mimetic discretisation *with* multipliers.
sols{2} = incompMimetic(xm, G, S1, fluid, 'Wells', Wm);

% Solve the system using TPFA *without* multipliers.
sols{3} = incompTPFA(xtp, G, T, fluid, 'Wells', Wtp);

% Solve the system using TPFA *with* multipliers.
sols{4} = incompTPFA(xtp, G, T .* m, fluid, 'Wells', Wtp);

%% Visualise the results
% We finally align the pressure results in a 2-by-2 matrix with the results
% unaffected by multipliers on the left and the multiplier results on the
% right.  The top row is discretised using the mimetic method while the
% bottom row shows the TPFA results.
clf
ttxt = {'Mimetic without multipliers', 'Mimetic with multipliers', ...
    'TPFA without multipliers', 'TPFA with multipliers'};
for i=1:4
    subplot(2,2,i)
    plotCellData(G, convertTo(sols{i}.pressure, barsa), ...
        'EdgeColor', 'k', 'EdgeAlpha', 0.1)
    view(-80,24), colorbar('horiz'), axis tight off
    set(gca, 'DataAspectRatio', daspctrat), zoom(1.9)
    title(ttxt{i})
end
for i=4:-1:1, subplot(2,2,i), set(gca,'Position',get(gca,'Position')-[0 .1 0 0]); end

%%
% We notice that the results are qualitatively the same between the
% discretisations (rows), while distinctly different between the models
% with and without multiplier effects.  This demonstrates that it is
% possible to include barrier modelling in the mimetic method by means of
% transmissibility multipliers, at least when the multipliers are between
% zero and one (i.e., $m \in (0,1)$), and that barrier effects can be quite
% dramatic.

%% Copyright notice
displayEndOfDemoMessage(mfilename)

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
