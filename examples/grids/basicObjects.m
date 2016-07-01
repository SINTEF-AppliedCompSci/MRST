%% Basic Objects in MRST Workflows
%
% The most fundamental object in MRST is the grid.  This dual purpose
% object both discretises a model's geometry and serves as a database to
% which static and dynamic properties and other objects are attached.
% Examples of such associate data are
%
% * Static petrophysical properties (e.g., permeability tensor or porosity)
% * Derived properties like transmissibility or pore-volume
% * Dynamic mass distributions (phase saturation and/or concentration of
%   chemical species or tracers)
% * Explicit model partitionings like PVT or relative permeability region
%   mappings
% * Source terms like wells or open boundaries
%
% This example explores a number of ways of constructing grids with degrees
% of complexity varying from simple Cartesian models with synthetic values
% for the porosity and permeability through to fully unstructured models in
% the industry-standard corner-point format.  We moreover introduce a
% collection of open datasets that are either hosted on the MRST website or
% published on other sites.

%% Activate Requisite Modules
%
% MRST comes equipped with a graphical user interface for interactively
% exploring static and dynamic properties attached to a grid.  We enable
% this functionality by activating the |mrst-gui| module that is provided
% in the MRST distribution.  Function <matlab:help('mrstModule')
% |mrstModule|> is the gateway to MRST's module system.  Modules are add-on
% features that can be enabled and disabled on demand.  For now we will not
% discuss this feature any further, but interested readers can find out
% more about it in the
% <http://www.sintef.no/projectweb/mrst/downloadable-resources/release-notes-for-mrst-2014a/#ModSysTutorial
% Release Notes for MRST 2014a>. We also activate the |deckformat| module,
% which contains functionality for reading and parsing data files given in
% the industry-standard ECLIPSE format.

mrstModule add mrst-gui deckformat

%% Simple Cartesian Model
%
% We start with a simple Cartesian model.  Function |cartGrid| creates
% shoe-box like models using a pair (for 2D) or triplet (for 3D) of
% integers representing the number of cells in each cardinal direction.
% Unless otherwise instructed, |cartGrid| will create cells of unit
% dimension (one metre in each direction).  Consequently, the total size of
% the model geometry depends on the number of cells.  This behaviour may be
% altered by passing the total size of the shoe-box as a second argument.
% The first grid is a 5-by-5-by-10 discretisation of the unit cube.  We can
% visualise the resulting geometry using function |plotGrid|.

G = cartGrid([5, 5, 10], [1, 1, 1].*meter);

plotGrid(G), view([-40, 32]), xlabel('x'), ylabel('y'), zlabel('Depth')

%% Populate Model
%
% Function |makeRock| creates a structure that holds the porosity and
% permeability fields, customarily referred to as |rock| throughout the
% MRST package.  Function |plotCellData| can then be employed to visualise
% the porosity field.

rock = makeRock(G, 100*milli*darcy, 0.3);

plotCellData(G, rock.poro)
view([-40, 32]), xlabel('x'), ylabel('y'), zlabel('Depth'), colorbar

%% Deactivate Selected Cells
%
% Function |gridLogicalIndices| will, provided the grid has an underlying
% Cartesian structure, compute (I,J,K) triplets either for all grid cells
% or a user-selected subset of those grid cells.  We can, for instance, use
% these values to pick out a subset of the cells and exclude those cells
% from the model using function |removeCells|.
%
% We'll pick out the 3-by-3-by-6 sub-cube in the centre of the model.  This
% cube can be identified by the (I,J) indices having values in the range
% 2:4 and the K indices being in the range 3:8.  Combining this with
% logical expressions means we can extract the subset in an economical way
% using MATLAB(R)'s built-in function |all|.

pick_ij = false([ 5, 1]);  pick_ij(2 : 4) = true;
pick_k  = false([10, 1]);  pick_k (3 : 8) = true;

ijk = gridLogicalIndices(G);
c   = all([pick_ij([ijk{[1, 2]}]), pick_k(ijk{3})], 2);
Gr  = removeCells(G, c);

clf
plot_cells = true([5, 1]);  plot_cells([1, 2]) = false;
plotGrid(G, 'FaceColor', 'none', 'EdgeAlpha', 0.6, ...
         'EdgeColor', [0.6, 0.4, 0.5], 'LineWidth', 1.2)
plotGrid(Gr, plot_cells(ijk{2}(~c)))
view([-48, 25]), xlabel('x'), ylabel('y'), zlabel('Depth')

material shiny, camlight(-60, -50)
for k = 1:3, camlight(-60, 50), end

%% Geostatistical Distribution of Petrophysical Properties
%
% MRST does not have advanced geostatistical functionality.  The two
% functions |logNormLayers| and |gaussianField| may be used to sample
% random fields in a spatial region with various degress of correlation,
% but we recommend that more complete geostatistical software be employed
% in a more realistic situations.  Here, we simply demonstrate that we can
% visualise spatially varying data.  Note that the formula that derives the
% permeability field from the porosity fields makes certain simplifying
% assumptions about the nature of the rock and its grains.  Function
% |compressRock| from the deckformat module is useful when extracting the
% petrophysical properties that pertain to a model's active cells only.

poro = reshape(gaussianField(G.cartDims, [0.2, 0.4]), [], 1);
perm = poro.^3 .* (1e-5)^2 ./ (0.81 * 72 * (1 - poro).^2);

rock   = makeRock(G, perm, poro);
rock_r = compressRock(rock, Gr.cells.indexMap);

clf, subplot(1,2,1)
plotCellData(G, log10(convertTo(rock.perm, milli*darcy)))
view([-48, 20]), xlabel('x'), ylabel('y'), zlabel('Depth')

subplot(1,2,2)
plotGrid(G, 'FaceColor', 'none', 'EdgeAlpha', 0.4, ...
         'EdgeColor', [0.6, 0.4, 0.5], 'LineWidth', 1.2)
plotCellData(Gr, log10(convertTo(rock_r.perm, milli*darcy)), ...
             plot_cells(ijk{2}(~c)))
view([-48, 20]), xlabel('x'), ylabel('y'), zlabel('Depth')

material shiny, camlight(-60, -50)
for k = 1:3, camlight(-60, 50), end

%% Unstructured Grids in Corner-Point Format
%
% MRST provides mature support for unstructured grids in the corner-point
% format as represented in the ECLIPSE reservoir simulator.  There are
% multiple example grids in the
% <matlab:what(fullfile(ROOTDIR,'gridprocessing','testgrids')) testgrids>
% directory.  Here, however, we will use one of the datasets known to MRST.
% Function |getAvailableDatasets| contains the list of known datasets.  The
% second return value, |present|, is a logical array that states whether or
% not the corresponding dataset is available on the local computer system.

[info, present] = getAvailableDatasets;

display([ { info(present).name }      ; ...
          { info(present).modelType } ] .')

%% Construct grid and properties from ECLIPSE-type input.
%
% Function |processGRDECL| constructs MRST grids from ECLIPSE-style input.
% Such data can be input using either function |readGRDECL| or function
% |readEclipseDeck|.  The latter is intended for complete simulation cases,
% i.e., ECLIPSE-style .DATA files, whereas the former reads loosely coupled
% grid data and does not form a state machine to handle temporally changing
% wells and related constraints.
%
% MRST's function |getDatasetPath| returns the full path of a particular
% dataset on the local computer system.  We take care to convert the input
% data to MRST's strict SI only unit conventions in the process.  Using
% strict SI units means that no solver or simulator needs to be aware of
% any particular unit convention.  Note moreover that in order to visualise
% the pore-volume field we need to compute the bulk (geometric) volume of
% all grid cells.  This calculation, along with cell and connection
% centroids, connection areas, and connection normals, is affected by
% MRST's built-in function |computeGeometry|.
grdecl = readGRDECL(fullfile(getDatasetPath('BedModel2'), ...
                             'BedModel2.grdecl'));

grdecl = convertInputUnits(grdecl, getUnitSystem('LAB'));

G    = computeGeometry(processGRDECL(grdecl));
rock = compressRock(grdecl2Rock(grdecl), G.cells.indexMap);
pvol = poreVolume(G, rock);

clf
plotCellData(G, convertTo(pvol, (centi*meter)^3), 'EdgeAlpha', 0.1), colorbar
view([-32, 22]), xlabel('x'), ylabel('y'), zlabel('Depth'), grid on

%%
% This model gives a high-resolution presentation of a 30×30×6 cm^3
% sedimentary bed that contains six different rock types. The model is
% represented in the form of a grid with 30×30×333 cells and is a good
% example of a corner-point grid having a large number of inactive cells
% and cells with degenerate geometry. We add the rock types to the |rock|
% structure and then launch the |plotToolbar| viewer to let you inspect the
% model in more detail. To better understand the model, you can, for
% instance, push the 'ijk' button to view the model in logical space rather
% than in physical space, or use the histogram picker to inspect the
% distribution of each individual rock type.

rock.rocktype = grdecl.SATNUM(G.cells.indexMap);
plotToolbar(G, rock,'field','rocktype');
view([-32, 22]), xlabel('x'), ylabel('y'), zlabel('Depth'), axis tight

%% Initialise Complete Simulation Case From ECLIPSE Data
%
% MRST is able to use complete ECLIPSE simulation cases as input and
% construct objects that are usable in an MRST simulation scenario.  This
% functionality is available in the deckformat module and centred around
% the functions |readEclipseDeck|, |convertDeckUnits|, |initEclipseGrid|,
% and |processWells|.  The wrapper function |initEclipseModel| calls these
% and other related functions in the appropriate order to create a fully
% initialised simulation model, complete with initial mass and pressure
% distributions.  We will visualise the initial phase distribution in the
% well-known Nineth SPE Comparative Solution Project, plotting the oil
% phase as green, the gas phase as red and the water phase as blue.

input = fullfile(getDatasetPath('SPE9'), 'BENCH_SPE9.DATA');

[G, rock, fluid, state, wells] = initEclipseModel(input);

name_wells = @(prefix, n) ...
   arrayfun(@(i) sprintf('$\\mathit{%s}_{%02d}$', prefix, i), ...
            1 : n, 'UniformOutput', false);

is_inj = [ wells.sign ] > 0;
wname  = name_wells('I', sum(is_inj));

[wells(is_inj).name] = wname{:};

is_prod = [ wells.sign ] < 0;
wname   = name_wells('P', sum(is_prod));

[wells(is_prod).name] = wname{:};

clf
plotCellData(G, state.s(:, [3, 2, 1])), axis tight, grid on
[htop, htext, hs, hline] = plotWell(G, wells);
view([55, 40]), xlabel('x'), ylabel('y'), zlabel('Depth')

%%
% The |ad-core| and |ad-blackoil| modules contain simulators and tools for
% simulating this type of models. In particular, the script
% |blackoilTutorialSPE9| discusses the various data elements that make up
% the SPE9 model in more detail, performs a full simulation, and analyzes
% the simulation results. To get a list of other tutorials that comes with
% these modules, type |mrstExamples ad-core ad-blackoil|.

%% Copyright Notice
%
% #COPYRIGHT_EXAMPLE#
