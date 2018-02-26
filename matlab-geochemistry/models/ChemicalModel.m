classdef ChemicalModel < PhysicalModel
    % ChemicalModel An object for the instantiation and solution of chemical
    % systems. 
    %
    % SYNOPSIS:
    %  chem = ChemicalModel(elementNames, speciesNames, reactions)
    %
    % DESCRIPTION:
    %   A class of PhysicalModel which can construct and solve
    %   arbitrary aqueous geochemical system including nonisothermal
    %   aqueous and surface speciation, redox chemistry, equilibrium 
    %   with solid and gas phases, and surface complexation.
    %
    % PARAMETERS:
    %   elementNames - A cell array of strings containing all
    %       elements to be considered in the chemical system. 
    %       Elements do not have to correspond to actual element names,
    %       but it is reccomened that they do. To specify a total element
    %       concentration as an input to the system append an asterisk, '*', to
    %       the element name.
    %            
    %           elementNames = {'H', 'O*'};
    %
    %   speciesNames - A cell array of strings of the chemical
    %       species to be considered in the chemical system. Each 
    %       entry in speciesNames must be a combination of
    %       entries found in elementNames, or surfaces. The species
    %       charge can be desginated with a +/- at the end of the
    %       name followed by the charge number (i.e. Ca+2, Cl-, >FeO-1/2).To 
    %       specify a species concentration as an input to the system append 
    %       an asterisk, '*', to the species name.
    %
    %           speciesNames = {'H+*', 'OH-', 'H2O'};
    %
    %   reactions - A cell array listing the chemical
    %       reactions of the system as strings, followed
    %       by the equilibrium constant of the reaction as a scalar
    %       quantity. Entries in reactions must be given as a
    %       string in the form 'reactants = products'. 
    %       Equilibrium constant must be given in SI units (mol/m^3).
    %
    %           reactions = {'H2O = H+ + OH-', 1e-14*mol/litre};
    %
    % KEYWORD ARGUMENT:
    %   'surfaces' - A cell structure containing the names
    %       of surface functional groups, followed by information
    %       regarding the specific surface. The first entry of
    %       surfaces is the name of the surface functional group.
    %       Surface functional group names must begin with the ">"
    %       symbol (i.e. '>FeO'). There are three different surface
    %       chemistry models implemented in ChemicalModel: the
    %       langmuir, the triple layer, and the
    %       constant capacitance models. The diffuse layer and basic
    %       stern models can be constructed by an adjustment of the
    %       capacitance parameters of the triple layer model. A chemical system
    %       can have any number of surfaces with different surface
    %       complexation models. The total concentration of surface sites will
    %       be calulated from the speciric surface area, site density and
    %       slurry density. 
    %
    %       The Langmuir model is the most simplistic case.
    %
    %           geometry = [1/(nano*meter)^2 50*meter^2/gram 1000*grams/litre];
    %
    %           surfaces = {'>FeO', {geometry, 'langmuir'}};
    %
    %           chem = ChemicalModel(elementNames, speciesNames, reactions,...
    %                       'surf', surfaces);
    %
    %       In this example the surface functional group is
    %    	'>FeO'. The second entry is a cell containing the
    %     	surface parameters. The variable geometry must
    %   	always be the first entry in this parameter cell.
    %     	The entries of geometry must correspond to the
    %      	surface site density, the specific surface area, and
    %      	the slurry density in that specific order. Values
    %     	must be given in SI units. The second entry of the
    %     	parameters cell must be the type of surface model to
    %     	used. In the case of 'langmuir' as above
    %     	no other entries are needed. 
    %
    %       The triple layer model requires additional information.
    %
    %           capacitances = [1 0.2]*Farad/meter^2;
    %
    %           surfaces = {'>FeO', {geometry, 'tlm', capacitance}
    %
    %           chem = ChemicalModel(elementNames, speciesNames, reactions,...
    %                       'surf', surfaces);
    %
    %       In this example geometry is defined as before, but now
    %       the '>FeO' is a triple layer model surface, as is
    %       indicated by 'tlm.' Additional parameters are needed in
    %       this case. The capacitance density of each of the
    %       Helmholtz layers must be given after the model type.
    %       The basic stern model can be simulated by increasing
    %       the capacitance of the second layer (i.e. [1 1e3]) and
    %       the diffuse layer model can be simulated by increasing
    %       the capacitance of both layers (i.e. [1e3 1e3]). If a
    %       surface species charge is distributed between multiple surfaces
    %       or contributes to a surface other than the mineral surface
    %       the charging behavior can be listed after capacitance values
    %
    %           surfaces = {'>FeO', {geometry, 'tlm', capacitance,...
    %                                           '>FeONa', [-1 +1 0]}
    %
    %       in the above example the species >FeONa contributes -1 charge to
    %       the surface and +1 charge to the inner Helmholtz plane. 
    %
    %       The constant capacitance model can be similarly defined.
    %
    %           capacitances = 1*Farad/meter^2;
    %
    %           surfaces = {'>FeO', {geometry, 'ccm', capacitance, '>FeO-'}
    %
    %           chem = ChemicalModel(elementNames, speciesNames, reactions,...
    %                       'surf', surfaces);
    %
    %       In this example geometry is defined as before, but now
    %       the '>FeO' is a constant capacitance surface, as is
    %       indicated by 'ccm.' The capacitance density of the
    %       Helmholtz layers must be given after the model type.
    %
    %       If surface funcational groups reside on the same surface they can be
    %       combined into a surface group, meaning all species will contribute
    %       to the charge conservation equations and share the same
    %       electrostatic properties of the group.
    %
    %           surfaces = {'>FeO', {FeOgeo, 'ccm', capacitance, '>FeO-'},...
    %                       '>Fe3O' {Fe3Ogeo, 'ccm', capacitance, '>FeO-'},...
    %                       'groups', {'iron', {'>FeO', '>Fe3O'});
    %
    %           chem = ChemicalModel(elementNames, speciesNames, reactions,...
    %                       'surf', surfaces);
    %
    %       Here the surface functional groups belong to the surface group
    %       'iron.' The name of the surface group can be any string. Multiple
    %       groups can be defined in this way. The electrostatic properties
    %       (capacitance) must be the same among surface functional groups that
    %       reside on the same surface group. 
    %
    %   'combinations' - a cell containing name value pairs of user defined
    %       linear combinations of species. This is how one would define 
    %       alkalinity or other such parameters. To specify a linear combination 
    %       concentration as an input to the system append an asterisk, '*', 
    %       to the combination's name.
    %
    %           combos = {'alk*', 'HCO3- + 2*CO3-2 + OH- - H+',...
    %                     'pot', 'PO4-2 + HPO4-2'};
    %
    % RETURNS:
    %   chem - the chemical model object containing tools used internally for
    %   the solution of the chemical system, and some other information which
    %   is useful for the user in plotting and data manipulation. Some useful
    %   fields include:
    %       chem.speciesNames - names of the species names in the system 
    %       chem.elementNames - names of the elements names in the system
    %       chem.reacionConstants - the equilibrium constants of the chemical
    %           reactions in the same order the reactions were defined
    %       chem.combinationNames - a cell array of strings of the linear
    %           combination names
    %       chem.inputs - a cell array of strings containing the inputs in the
    %           order which they are to be provided to chem.initState
    %       chem.solidNames - name of the solid phases in the system
    %       chem.gasNames - names of the gas phases in the system
    %
    %   editable fields include:
    %       chem.plotIterations - toggle plotting of concentration values at
    %           each iteration, true or [false]
    %       chem.nonLinearMinIterations - minimum number of iterations for the
    %           nonlinear solver, [1]
    %       chem.nonLinearMaxIterations - maximum number of iterations for the nonlinear solver, [25] 
    %       chem.nonLinearTolerance - tolerance of the residual of the
    %           nonlinear system, [1e-12]
    %       chem.linearTolerance - tolerance of the residual of the linear
    %           system, for backslash, [1e-8]
    %       chem.linearMaxIterations - maximum number of iterations for the
    %           linear solver, [25]
    %       chem.chargeBalanceTolerance - tolerance of charge balance equation
    %           as fraction of total ion concentration, [0.05]
    %        
    %   
    % EXAMPLES:
    %   specifying inputs - MATLAB geoChemistry leverages the automatic
    %       differentiation capabilities of MRST. This means that any variable
    %       can be an unknown or an input. This makes the tool incredibly
    %       flexible and especially well suited for the utilization of 
    %       labortaory data. Inputs are chosen when the chemical system is
    %       initiated with ChemicalModel by appending an asterisk (*) in the
    %       component name. If an element is chosen as an input, the user
    %       provides the total concentration of that element as an input.
    %       Aqueous and surface species as well as combination components can
    %       can all be chosen as inputs. 
    %
    %           elementNames = {'H', 'O'};
    %
    %           speciesNames = {'H+*', 'OH-', 'H2O*'};
    %
    %           reactions = {'H2O = H+ + OH-', 1e-14*mol/litre};
    %
    %           chem = ChemicalModel(elementNames, speciesNames, reactions);
    %
    %       In the above example the aqueous concentration of H+ and H2O has
    %       been specified as an input. 
    %
    %       The system is solved by passing the input values to ChemicalModel.initState
    %       The order of the inputs can be found by 
    %
    %           chem.inputs
    %
    %       initState accepts arrays or vectors, meaning that variable sweeps
    %       can be easily done
    %
    %           n = 100;
    %           H = logspace(-2, -10, n)'
    %           H2O = ones(n,1);
    %           inputs = [H H2O]*mol/litre;
    %
    %           state = chem.initState(inputs)
    %
    % SEE ALSO:
    %   'initState', 'printChemicalSystem'

    %{
    Copyright 2009-2017 SINTEF Digital, Mathematics & Cybernetics.

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

 
    properties
        elementNames      % Names of the master components as first master component.
        logElementNames   % Names of the log master components (= 'log' in
                             % front of elementNames).

        MCind                % index vector of each master component
        nMC                  % number of master components

        speciesNames            % Names of the components.
        logSpeciesNames         % Names of the log components (= 'log' in
                             % front of speciesNames).

        compositionMatrix    % Composition for each components, in term of the
                             % master components.
        Cind                 % index vector of species
        nC                   % Number of components
        
        maxMatrices          % cells of matrices, one per components. Used to
                             % compute the upper bound when updating a component.

        chargeVector         % coefficients of charged species.
        reactionMatrix       % Coefficients for the reaction.
        reactionConstants    % Contains the K_i.
        logReactionConstants % Contains the log K_i.


        nR                   % Number of reactions
        rxns                 % string of reactions


        surfChargeMatrix     % matrix of surface charge contributions
        surfMaster           % index vector of master components that are surfaces

        surfFlag             % are there surface species 1 or 0
        surfaces             % surface groups for electrostatic calculations
        
        nSurf                % number of surface components
        nP                   % number of surface potential terms
        
        compositionReactionModel % model for computing the mas conservation and reaction equations together
        chargeBalanceModel   % model for computing the solution when charge balance is required
        
        chemicalInputModel   % Model to solve chemical equilibrium state from user
                              % input. see initState member function. variable
                              % is initialized in initSecondaryComponents.

        surfInfo            % information regarding the surfaces
        activityNames   % name of activity of species

        surfaceChargeNames         % name of surface charge variables
        surfacePotentialNames         % name of potentials
        
        surfaceActivityCoefficientNames % surface potential multipliers
        logSurfaceActivityCoefficientNames %  log surface potential multipliers
        
        allCharge           % vector for ensuring each reaciton is charge balanced
        
        inputs              % names of inputs the user must provide
        
        combinationNames    % names of linear combinations
        combinationMatrix   % combination for each linear combination in terms of species
        nLC                 % number of linear combinations
        
        aqueousConcentrationNames % name of total aqueous concentration variables
        surfaceConcentrationNames % name of total surface concentration variables
        
        solidNames          % name of species in the solid phase
        gasNames            % name of species in the gas phase

        logSolidNames          % log name of species in the solid phase
        logGasNames            % log name of species in the gas phase
        
        nG                  % number of gas components
        nS                  % number of solid components
        
        allContributionMatrix % matrix to check if all reactions are balanced
        
        allComponentNames        % name of all components, in all phases
        phaseInd            % index of phase names in component names
        gasInd              % index of gas phase in comp names
        solidInd            % indexof solid phase in comp names
        
        allReactionMatrix   % matrix of contribution of all components in all phases in a reaction
        partialPressureNames % names of partial pressures of gasses
        solidDensityNames % names of solid densities
        
        allCombinationMatrix   % combination matrix for printing the chemical system

        gasContributionMatrix       % composition matrix of gas phase
        solidContributionMatrix     % composition matrix of solid phase
        
        gasReactionMatrix   % reaction matrix of gas phase
        solidReactionMatrix % reaction matrix of solid phase
                
        surfacePotentialMatrix  % matrix for surface potential contributions to reactions
        
        plotIterations            % toggle plotting of iteration information true or false
        nonLinearMinIterations    % minimum number of iterations for the nonlinear solver
        nonLinearMaxIterations    % maximum number of iterations for the nonlinear solver  
        nonLinearTolerance        % tolerance of the residual of the nonlinear system
        linearTolerance           % tolerance of the residual of the linear system, for backslash
        linearMaxIterations       % maximum number of iterations for the linear solver
        
        CVC                     % charge variation component
    end


    methods

        function model = ChemicalModel(varargin)
 
            
    
            model = model@PhysicalModel([]);
            if nargin >= 3
                
                % parse inputs to chemical model
                masterComponentNames = varargin{1};
                componentNames       = varargin{2};
                reactions       = varargin{3};
                
            	p = inputParser;
            
                valFun = @(x) iscell(x);
                p.addParameter('surfaces', '', valFun);                
                p.addParameter('combinations', '', valFun);

                p.parse(varargin{4:end})
            
                % add surface functional groups to master component list
                [model, toPass, masterComponentNames] = parseSurfaceComponents(model, masterComponentNames, p);
                
                
                % initiate master components
                model = initMasterComponents(model, masterComponentNames);
                % initiate species
                model = initSecondaryComponents(model,  componentNames);
                % initiate electrostatic model and components
                model = initElectrostaticModel(model, toPass );
                % initiate linear combinations
                model = initLinearCombinations(model, p);

                
                % initiate chemical reactions
                model = initChemicalReactions(model, reactions);
                
                % initiate max value matrix for variables bounding
                model = setupMaxMatrices(model);
                
                % instantiate the chemical input solver model and parse
                % inputs
                props = properties(model);
                for i = 1 : numel(props);
                    if ~strcmpi('compositionReactionModel', props{i})
                        model.chemicalInputModel.(props{i}) = model.(props{i});
                    end
                end
                model.inputs = model.chemicalInputModel.inputNames;
                
                
                if model.surfFlag
                    inInd = zeros(1,numel(model.inputs));
                    for i = 1 : numel(model.surfInfo.master)
                        inInd = inInd + strcmpi(model.chemicalInputModel.inputNames, model.surfInfo.master{i});
                    end
                    model.inputs = model.chemicalInputModel.inputNames(~inInd);
                end
                                
                
            elseif nargin == 0
                % Used when creating instance of
                % chemicalInputModel. Nothing is initiated at this stage.
            else
                error(['The constructor of ChemicalModel is called without the ', ...
                       'appropriate number of input parameters']);
            end

            model.plotIterations            = false;
            model.nonLinearMinIterations    = 1; 
            model.nonLinearMaxIterations    = 25;  
            model.nonLinearTolerance        = 1e-12;
            model.linearMaxIterations       = 25;
            model.linearTolerance           = 1e-8;
            model.CVC = '';
        end


        %%
        function [model, toPass, masterComponentNames] = parseSurfaceComponents(model, masterComponentNames, p)
            
            if isempty(p.Results.surfaces)
                toPass = '';
                return
            end
            
            assert(rem(numel(p.Results.surfaces),2)==0, 'A key/value pair must be provided for surfaces, yet the number of inputs is odd.');

            surfNames = p.Results.surfaces(1:2:end);
            ind = cellfun(@(x) ~isempty(x), regexpi(surfNames,'g[roups]'));
            surfNames = surfNames(~ind);

            %%% check which surface masters are in the defined groups and
            %%% then remove them from groups as above
            if sum(ind) > 0
                tmp = find(ind);
                defGroups = p.Results.surfaces{2*tmp};
                defGroupNames = defGroups(1:2:end);
                defGroupMast = defGroups(2:2:end);

                toPass = p.Results.surfaces;
                toPass(2*tmp) =[];
                toPass(2*tmp-1) = [];

                % make sure group names are unique
                for i = 1 : numel(defGroupNames)
                    ind = zeros(size(defGroupNames));
                    ind(i) = 1;
                    test = true;
                    if numel(defGroupNames) > 1
                        test = ~all(strcmpi(defGroupNames{i}, defGroupNames(~ind)));
                    end
                    assert(test, 'Surface group names must be unique.');
                end

                % make sure surface functional groups do not appear in more than one
                % user defined group
                spnames = {};
                for i = 1 : numel(defGroupMast);
                    spnames = [spnames defGroupMast{i}];
                end

                for i = 1 : numel(spnames)
                    ind = zeros(size(spnames));
                    ind(i) = 1;

                    test = all(~strcmpi(spnames{i}, spnames(~ind)));

                    validatestring(spnames{i}, surfNames, 'ChemicalModel', 'entries to share');
                    assert(test, ['The surface species ' spnames{i} ' is repeated or appears in more than one group.']);
                end

                % add functional groups not listed in groups to
                % groupNames may need to change cell indexing for when there
                % arent the same number of surface master components in  a group
                groupNames = defGroupNames;
                masterNames = defGroupMast;
                for i = 1 : numel(surfNames);
                    if ~any(strcmpi(surfNames{i}, spnames))
                        groupNames{end+1} = surfNames{i};
                        masterNames{end+1} = surfNames(i);
                    end
                end
            else
                for i = 1:numel(surfNames)
                    groupNames = surfNames;
                    masterNames{i} = surfNames(i);
                end
            end

            model.surfaces.groupNames = groupNames;
            model.surfaces.masterNames = masterNames;

            surfNames = cellfun(@(name) [name '*'], ...
                                       surfNames, ...
                                       'uniformoutput', false);

            masterComponentNames = horzcat(masterComponentNames, surfNames);
            
            if ~exist('toPass', 'var')
                toPass = p.Results.surfaces;
            end
                
        end
                          
        %%
        function [problem, state] = getEquations(model, state0, state, dt, drivingForces, varargin)

            [logcomps, logmasterComps] = prepStateForEquations(model, state);

            [eqs, names, types] = equationsChemicalLog(logcomps, logmasterComps, ...
                                                       model);

            primaryVariables = model.logSpeciesNames;
            problem = LinearizedProblem(eqs, types, names, primaryVariables, state, dt);

        end



        %%
        function [fn, index] = getVariableField(model, name)

            varfound = false;

            while ~varfound


                if strcmpi(name, 'partialPressures')
                    varfound = true;
                    fn = 'partialPressures';
                    index = ':';
                    break
                end

                ind = strcmpi(name, model.partialPressureNames);
                if any(ind)
                    varfound = true;
                    fn = 'partialPressures';
                    index = find(ind);
                    break
                end

                if strcmpi(name, 'fluidVolumeFraction')
                    varfound = true;
                    fn = 'fluidVolumeFraction';
                    index = ':';
                    break
                end
                
                if strcmpi(name, 'logFluidVolumeFraction')
                    varfound = true;
                    fn = 'logFluidVolumeFraction';
                    index = ':';
                    break
                end

                 if strcmpi(name, 'matrixVolumeFraction')
                    varfound = true;
                    fn = 'matrixVolumeFraction';
                    index = ':';
                    break
                end
                
                if strcmpi(name, 'logMatrixVolumeFraction')
                    varfound = true;
                    fn = 'logMatrixVolumeFraction';
                    index = ':';
                    break
                end
                
                if strcmpi(name, 'species')
                    varfound = true;
                    fn = 'species';
                    index = ':';
                    break
                end

                if strcmpi(name, 'chargeBalance')
                    varfound = true;
                    fn = 'chargeBalance';
                    index = ':';
                    break
                end
                
                if strcmpi(name, 'elements')
                    varfound = true;
                    fn = 'elements';
                    index = ':';
                    break
                end

                if strcmpi(name, 'logSpecies')
                    varfound = true;
                    fn = 'logSpecies';
                    index = ':';
                    break
                end
                
                if strcmpi(name, 'activities')
                    varfound = true;
                    fn = 'activities';
                    index = ':';
                    break
                end

                if strcmpi(name, 'logElements')
                    varfound = true;
                    fn = 'logElements';
                    index = ':';
                    break
                end

                if strcmpi(name, 'combinationComponents')
                    varfound = true;
                    fn = 'combinationComponents';
                    index = ':';
                    break
                end
                 
                ind = strcmpi(name, model.combinationNames);
                if any(ind)
                    varfound = true;
                    fn = 'combinationComponents';
                    index = find(ind);
                    break
                end
              
                
                ind = strcmpi(name, model.speciesNames);
                if any(ind)
                    varfound = true;
                    fn = 'species';
                    index = find(ind);
                    break
                end
                
                if strcmpi(name, 'aqueousConcentrations')
                    varfound = true;
                    fn = 'aqueousConcentrations';
                    index = ':';
                    break
                end

                ind = strcmpi(name, model.aqueousConcentrationNames);
                if any(ind)
                    varfound = true;
                    fn = 'aqueousConcentrations';
                    index = find(ind);
                    break
                end
                
                if strcmpi(name, 'surfaceConcentrations')
                    varfound = true;
                    fn = 'surfaceConcentrations';
                    index = ':';
                    break
                end

                ind = strcmpi(name, model.surfaceConcentrationNames);
                if any(ind)
                    varfound = true;
                    fn = 'surfaceConcentrations';
                    index = find(ind);
                    break
                end

               if strcmpi(name, 'partialPressures')
                    varfound = true;
                    fn = 'partialPressures';
                    index = ':';
                    break
                end

                ind = strcmpi(name, model.gasNames);
                if any(ind)
                    varfound = true;
                    fn = 'partialPressures';
                    index = find(ind);
                    break
                end

                if strcmpi(name, 'saturationIndicies')
                    varfound = true;
                    fn = 'saturationIndicies';
                    index = ':';
                    break
                end

                ind = strcmpi(name, model.solidNames);
                if any(ind)
                    varfound = true;
                    fn = 'saturationIndicies';
                    index = find(ind);
                    break
                end
                
                ind = strcmpi(name, model.solidDensityNames);
                if any(ind)
                    varfound = true;
                    fn = 'solidDensities';
                    index = find(ind);
                    break
                end
                
                 if strcmpi(name, 'solidDensities')
                    varfound = true;
                    fn = 'solidDensities';
                    index = ':';
                    break
                end
                
                if strcmpi(name, 'logPartialPressures')
                    varfound = true;
                    fn = 'logPartialPressures';
                    index = ':';
                    break
                end

                ind = strcmpi(name, model.logGasNames);
                if any(ind)
                    varfound = true;
                    fn = 'logPartialPressures';
                    index = find(ind);
                    break
                end

                if strcmpi(name, 'logSaturationIndicies')
                    varfound = true;
                    fn = 'logSaturationIndicies';
                    index = ':';
                    break
                end

                ind = strcmpi(name, model.logSolidNames);
                if any(ind)
                    varfound = true;
                    fn = 'logSaturationIndicies';
                    index = find(ind);
                    break
                end
                
                ind = strcmpi(name, model.activityNames);
                if any(ind)
                    varfound = true;
                    fn = 'activities';
                    index = find(ind);
                    break
                end
                
                ind = strcmpi(name, model.logSpeciesNames);
                if any(ind)
                    varfound = true;
                    fn = 'logSpecies';
                    index = find(ind);
                    break
                end
                
                ind = strcmpi(name, model.elementNames);
                if any(ind)
                    varfound = true;
                    fn = 'elements';
                    index = find(ind);
                    break
                end

                ind = strcmpi(name, model.logElementNames);
                if any(ind)
                    varfound = true;
                    fn = 'logElements';
                    index = find(ind);
                    break
                end

                ind = strcmpi(name, 'temperature');
                if any(ind)
                    varfound = true;
                    fn = 'temperature';
                    index = find(ind);
                    break
                end
              
                
                ind = strcmpi(name, model.surfaceChargeNames);
                if any(ind)
                    varfound = true;
                    fn = 'surfaceCharges';
                    index = find(ind);
                    break
                end
                
                ind = strcmpi(name, model.surfacePotentialNames);
                if any(ind)
                    varfound = true;
                    fn = 'surfacePotentials';
                    index = find(ind);
                    break
                end
                
                ind = strcmpi(name, model.surfaceActivityCoefficientNames);
                if any(ind)
                    varfound = true;
                    fn = 'surfaceActivityCoefficients';
                    index = find(ind);
                    break
                end
                
                if strcmpi(name, 'surfaceActivityCoefficients')
                    varfound = true;
                    fn = 'surfaceActivityCoefficients';
                    index = ':';
                    break
                end
                
                ind = strcmpi(name, model.logSurfaceActivityCoefficientNames);
                if any(ind)
                    varfound = true;
                    fn = 'logSurfaceActivityCoefficients';
                    index = find(ind);
                    break
                end
                
                if strcmpi(name, 'logSurfaceActivityCoefficients')
                    varfound = true;
                    fn = 'logSurfaceActivityCoefficients';
                    index = ':';
                    break
                end
                
                [fn, index] = getVariableField@PhysicalModel(model, name);
                varfound = true;
            end

        end

        %%
        function fds = getAllVarsNames(model)
            fds = {'species', 'elements', 'logSpecies', ...
                   'logElements', model.speciesNames{:}, ...
                   model.elementNames{:}, model.logSpeciesNames{:}, ...
                   model.logElementNames{:}, 'logFluidVolumeFraction', model.logSurfaceActivityCoefficientNames{:}, model.logGasNames{:}, model.logSolidNames{:}};
        end

        %%
        function [state, report, model] = initState(model, userInput, varargin)
        % initState solves the chemical system using the specified input
        %
        % SYNOPSIS:
        %  [state, report, model] = initState(model, userInput, varargin)
        %
        % DESCRIPTION:
        %   A function used to solve the chemical system generated by
        %   ChemicalModel according to constraints given in userInput. 
        %
        % REQUIRED PARAMETERS:
        %   userInput - A matrix containing the values of elements
        %               or species specified as inputs during the generation of the
        %               chemical system using ChemicalModel. userInput must have as
        %               many columns as there are specified knowns, which is also equal
        %               to the number of entries in elements. 
        %           
        %
        % KEYWORD ARGUMENTS:
        %   'temperature' - 'key'/value pair representing the
        %                   temperature of the chemical system. Can be a scalar, or a
        %                   vector with the same length as userInput. If temperature is
        %                   not provided 298 K is used. 
        %          
        %    	state = model.initState(userInput, 'temp', 305)
        %
        %   'chargeBalance' - Toggles the enforcment of strict
        %                     charge balance. A master component must be specified as a 
        %                     'key'/value pair. If charge balance is not satisfied with the
        %                     given constraint the solver will adjust the total concentration
        %                     of the specified element to enforce charge balance. The 
        %                     specified element must be one of the elements provided as an
        %                     input. By defualt charge balance is not strictly enforced. 
        %                     This make the chemical solver more robust. 
        %
        %    	state = model.initState(userInput, 'chargeBalance', 'Cl')
        %
        %
        % RETURNS:
        %   state - A structure containing the solution to the
        %           chemical system. This includes state.elements and
        %           state.species. state.elements contains the total
        %           element concentrations and surface functional group
        %           concentrations. state.species contains the species
        %           concentrations. Values can be retreived from state using the
        %           model.getProp() function.
        %
        %   	OH = chem.getProp(state, 'OH-');
        %
        %
        %   report - The nonlinear solver report for the chemical system. 
        %
        %   model - The updated chemical model object.
        %
        % SEE ALSO:
        %   'ChemicalModel', 'getProps'
        
        %{
        Copyright 2009-2017 SINTEF Digital, Mathematics & Cybernetics.

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
            % parse inputs to initState
            p = inputParser;
            
            valFun = @(x) ischar(x);
            p.addParameter('chargeBalance', 'nochargebalance', valFun);
            p.addParameter('temperature', 298, @isnumeric);

            
            p.parse(varargin{:})
            
            if strcmpi('nochargebalance',p.Results.chargeBalance)
               
            elseif any(ismember({'e','e-'},p.Results.chargeBalance))
                error('Electron concentration, e or e-, can not be used for charge balance.')
                
            elseif sum(ismember(model.elementNames,p.Results.chargeBalance))==0
                warning('Using any quantitiy other than a total element concentration as a charge balance can yield unexpected results.')
                
            end
            
                
            nI = size(userInput,1);
            
            % grab temperature
            nTemp = size(p.Results.temperature, 1);
            if nTemp == 1
                state.temperature = repmat(p.Results.temperature, nI, 1);
            else
                assert(nTemp == nI, 'The number of cells in state.temperature do not correspond to the size of userInput. Repmat the input vector if you want to do a temperature sweep.');
                state.temperature = p.Results.temperature(:);
            end

            
            givenTest = any(strcmpi(p.Results.chargeBalance, horzcat(model.chemicalInputModel.inputNames,'nochargebalance')));
            if ~givenTest
                warning 'Using a total element concentration whos values are not given (marked with "*") for charge balance may result in unexpected results.';
            end
            
            chargeBalance = ~strcmpi(p.Results.chargeBalance,'nochargebalance');

            model.chemicalInputModel = model.chemicalInputModel.validateModel();
            
            props = properties(model);
            for i = 1 : numel(props);
                model.chemicalInputModel.(props{i}) = model.(props{i});
            end
            
            
            if isempty(model.compositionReactionModel)
                model.compositionReactionModel = compositionReactionModel();
            end
            
            props = properties(model.chemicalInputModel);
            for i = 1 : numel(props);
                model.compositionReactionModel.(props{i}) = model.chemicalInputModel.(props{i});
            end
                
            model.compositionReactionModel = model.compositionReactionModel.validateModel();

            cheminput = model.chemicalInputModel;
            inputNames = cheminput.inputNames;

            inSize = size(userInput);
            if ~isempty(model.surfInfo);
                k = numel(model.surfInfo.master);
            else
                k = 0;
            end
            
            assert(inSize(2) == model.nMC-k, ['For the specified chemical system the input constraint must have ' ...
                                num2str(model.nMC-k) ' columns.']);
            
            state.elements      = mol/litre*ones(size(userInput,1), model.nMC);
            state.species            = mol/litre*ones(size(userInput,1), model.nC);
            
            state.combinationComponents = ones(size(userInput,1), model.nLC);
            
            state.partialPressures         = 0.1*ones(size(userInput,1), model.nG);
            state.logPartialPressures      = log(state.partialPressures);
            
            state.saturationIndicies       = 0.1*ones(size(userInput,1), model.nS);            
            state.logSaturationIndicies    = log(state.saturationIndicies);
    
 
            % put the value of the surface functional groups into the input
            % vector
            call = 0;
            for i = 1 : model.nMC
                
                if regexpi(inputNames{i}, '>')
                    sInd = strcmpi(inputNames{i}, model.surfInfo.master);
                    d = model.surfInfo.d{sInd};
                    s = model.surfInfo.s{sInd};
                    a = model.surfInfo.a{sInd};
                    in{i} = d.*s.*a;
                    call = call + 1;
                else
                    in{i} = userInput(:,i-call);
                end
                
            end
            
            % make sure input is the correct size
            givenSize = cellfun(@(x) size(x,1), in);
            mSize = max(givenSize);
            assert(any(givenSize == mSize | givenSize == 1), 'Input size must be consistent between input and surface parameters');
           
            % set their values inside of state
            for i = 1 : model.nMC
            	state = model.setProp(state, inputNames{i}, in{i}.*ones(mSize,1));
            end
            
            state = model.syncLog(state);
            
            % create initial guess
            fprintf('Computing initial guess...\n')
            [state, ~, report] = model.compositionReactionModel.solveChemicalState(state);
            
            % solve chemical system
                
            if chargeBalance
                fprintf('Solving the chemical system with strict charge balance...\n');
                state0 = state;
                if isempty(model.chargeBalanceModel)
                    model.chargeBalanceModel = chargeBalanceModel();
                end
                for i = 1 : numel(props);
                    model.chargeBalanceModel.(props{i}) = model.chemicalInputModel.(props{i});
                end
                model.chargeBalanceModel.CVC = p.Results.chargeBalance;
                [state, ~, report] = model.chargeBalanceModel.solveChemicalState(state);
                
                if ~report.Converged
                    state = state0;
                    warning('Charge balance not acheived with given inputs, consider choosing a different component for charge compensation. Attempting solution without charge balance.');
                end
            end
            if ~report.Converged || ~chargeBalance
                fprintf('Solving chemical system...\n')
                [state, ~, report] = model.chemicalInputModel.solveChemicalState(state);
                if ~report.Converged
                    warning('Solver did not converge, use results with caution.');
                end
            end
            
            
        end

        %%
        function model = initMasterComponents(model, names)
            %initMasterSpecies parses input of master variable names in geochemical system
            %for later construction of secondary species.
          

            %remove white space
            names = regexprep(names,'[\s]','');

            % make sure input is a cell
            assert(iscell(names),'input must be a cell array of strings');
            
            test = cellfun(@(x) ~isempty(x), regexpi(names, '(psi)(sig)'));
            assert(sum(test) == 0, 'The use of "sig" and "psi" are reserved for use within ChemicalModel, please remove them from the element cell');
            % find mass input species and remove the asterik from the name
            InInd = cellfun(@(x) ~isempty(x), regexp(names, '*'));
            names = regexprep(names, '*','');

            if isempty(model.chemicalInputModel)
                model.chemicalInputModel = ChemicalInputModel();
            end
            model.chemicalInputModel.inputNames = ...
                horzcat(model.chemicalInputModel.inputNames, names(InInd));

            % look for surface species
            model.surfMaster = double(cellfun(@(x) ~isempty(x), regexp(names, '^>'))');
            model.surfFlag = sum(model.surfMaster) > 0;
            model.nSurf = sum(model.surfMaster);

            % check that master component names are not identical
            nM = numel(names);

            for i = 1 : nM
                assert(~(sum(strcmp(names{i},names)) > 1), ['Element ' ...
                                    'names must be unique. Surface functional groups should not be included in surfaces, not elements.']);
            end

            % create the index matricies for construction of secondary species
            tmpvec = zeros(nM, 1);
            for i = 1 : nM
                model.MCind{i} = tmpvec;
                model.MCind{i}(i) = 1;
            end

            % make a list of the component names
            model.elementNames = cell(1,nM);
            [model.elementNames{:}] = deal(names{:});

            model.logElementNames = cellfun(@(name) ['log', name], ...
                                               model.elementNames, ...
                                               'uniformoutput', false);

            % write number of master components
            model.nMC = nM;

        end

        %%
        function [model] = initSecondaryComponents(model, secondaryComponents )
            %initSecondarySpecies computes species as linear combinations
            %of master components
         

            assert(~isempty(model.elementNames),'The function initSecondaryComponents can only be called after initMasterComponents.');

            % make sure input is a name value pair of the correct size
            assert(iscell(secondaryComponents),'input must be a cell array of strings');

            names = regexprep(secondaryComponents,'[\s]','','ignorecase');

            test = cellfun(@(x) ~isempty(x), regexpi(secondaryComponents, '(psi)(sig)(_aq)(_surf)'));
            assert(sum(test) == 0, 'The use of "sig", "psi", "_aq", and "_surf" are reserved for use within ChemicalModel, please remove them from the species names.');
            
            % find input components and remove the asterik from the name
            InInd = cellfun(@(x) ~isempty(x), regexp(names, '*'));
            names = regexprep(names, '*','');
            
            % find phase components
            solidInd = cellfun(@(x) ~isempty(x), regexp(names, '(s)'));
            gasInd = cellfun(@(x) ~isempty(x), regexp(names, '(g)'));
            
            phaseInd = logical(solidInd + gasInd);
            model.phaseInd = phaseInd;
            model.gasInd = gasInd;
            model.solidInd = solidInd;
            
            if isempty(model.chemicalInputModel)
                model.chemicalInputModel = ChemicalInputModel();
            end
            model.chemicalInputModel.inputNames = ...
                horzcat(model.chemicalInputModel.inputNames, names(InInd));


            nI = numel(model.chemicalInputModel.inputNames);
            namesO = names;
              
            % unwrap master component vectors
            nS = numel(names);
            nM = model.nMC;

            % make sure each secondary component name is distinct from the
            % others and master components
            for i = 1 : nM
                assert(~(sum(strcmpi(names{i},names)) > 1), 'Species names must be unique.');
                assert(~(sum(strcmpi(names{i}, model.elementNames)) > 0), ['The chemical species ' '"' names{i} '" is a repeat of an element or surface functional group.']);
            end

            names = strrep(names, '(s)','');
            names = strrep(names, '(g)','');
            
            % length of master component string
            strlen = cellfun(@length,model.elementNames);

            % sort strings by length to avoid partial overwrites
            [~, I] = sort(strlen, 2, 'descend');

            % capture the charge
            tmpvec = zeros(1, nS);
            tmp = regexp(names, '+(\d*([./]\d*)?)*', 'tokens'); 
            for i = 1 : numel(tmp)
                if ~isempty(tmp{i})
                    if strcmp('', tmp{i}{1})
                        tmpvec(i) = 1;
                    else
                        try
                            tmpvec(i) = str2num(tmp{i}{1}{1});
                        catch
                            error(['Unrecognized charge for species ' names{i} '.']);
                        end
                    end
                end
            end
           tmp = regexp(names, '-(\d*([./]\d*)?)*', 'tokens'); 
            for i = 1 : numel(tmp)
                if ~isempty(tmp{i})
                    if strcmp('', tmp{i}{1})
                        tmpvec(i) = -1;
                    else
                        try
                            tmpvec(i) = -str2num(tmp{i}{1}{1});
                        catch
                            error(['Unrecognized charge for species ' names{i} '.']);
                        end
                    end
                end
            end

            
            model.allCharge = tmpvec;
            model.chargeVector = tmpvec;
            surfInd = cellfun(@(x) ~isempty(x), regexpi(namesO,'>'));
            model.chargeVector(surfInd) = 0;

            % check that gas and solid components do not have a charge
            assert(all(model.allCharge(solidInd)==0),'Species in the solid phase, (s), can not have a charge.');
            assert(all(model.allCharge(gasInd)==0),'Species in the gas phase, (g), can not have a charge.');

            model.chargeVector(:,phaseInd) = [];
            
            % replace master component names with index vectors
            for i = 1 : nS;
                names{i} = regexprep(names{i}, '[-+](\d*([./]\d*)?)*','','ignorecase');
                for j = 1 : nM
                    ind = I(j);
                    [match,nomatch] = regexpi(names{i}, '(model\.MCind{\d*})','match','split');
                    nomatch = strrep(lower(nomatch), lower(model.elementNames{ind}), ['model.MCind{' num2str(ind) '}']);
                    names{i} = strjoin(nomatch,match);
                end
                names{i} = regexprep(names{i}, '}(\d+)m' , '}*$1+m','ignorecase');
                names{i} = regexprep(names{i}, '}(\d+)$' , '}*$1','ignorecase');
                names{i} = regexprep(names{i}, '}m' , '}+m' , 'ignorecase');
                names{i} = regexprep(names{i}, '}(\d+)+' , '}*$1','ignorecase');
                names{i} = regexprep(names{i}, '}(' , '}+(','ignorecase');
                
                names{i} = regexprep(names{i}, ')(\d+)m' , ')*$1+m','ignorecase');
                names{i} = regexprep(names{i}, ')(\d+)$' , ')*$1','ignorecase');
                names{i} = regexprep(names{i}, ')m' ,       ')+m' , 'ignorecase');
                names{i} = regexprep(names{i}, ')(\d+)+' , ')*$1','ignorecase');
                names{i} = regexprep(names{i}, ')(' , ')+(','ignorecase');
                names{i} = regexprep(names{i}, '(\d+)(' , '$1+(','ignorecase');
            end
            

            tmpvec = zeros(nS, 1);

            for i = 1 : nS
                model.Cind{i} = tmpvec;
                model.Cind{i}(i) = 1;
            end
            


            % store cell array of secondary component names
            model.allComponentNames = cell(1,nS);
            [model.allComponentNames{:}] = deal(namesO{:});            
            
            nS = nS - sum(phaseInd);
             
             % store cell array of secondary component names
            model.speciesNames = cell(1,nS);
            [model.speciesNames{:}] = deal(namesO{~phaseInd});

            % store cell array of gas component names
            model.gasNames = cell(1,sum(gasInd));
            [model.gasNames{:}] = deal(namesO{gasInd});
            
            % store cell array of solid component names
            model.solidNames = cell(1,sum(solidInd));
            [model.solidNames{:}] = deal(namesO{solidInd});
            
            for i = 1 : numel(model.gasNames)
                model.partialPressureNames{i} = ['p' model.gasNames{i}];
            end
            
            for i = 1 : numel(model.solidNames)
                model.solidDensityNames{i} = ['d' model.solidNames{i}];
            end
            
            % create the mass constraint matrix
            try
                indsum = cell(1, numel(names));
                for i = 1 : numel(names)
                    indsum{i} = eval([names{i} ,';']);
                    assert(isnumeric(indsum{i}),['The chemical species "' namesO{i} '" appears to contain a string combination that does not correspond to an element or surface functional group. Make sure surface species contain ">" at the begining and are consistent with the surface name in surfaces.']);
                end
            catch
               error(['The chemical species "' namesO{i} '" appears to contain a string combination that does not correspond to an element or surface functional group. Make sure surface species contain ">" at the begining and are consistent with the surface name in surfaces.']);
            end
            model.compositionMatrix = horzcat(indsum{:});

            % remove contribution of solid to composition matrix
            model.allContributionMatrix = model.compositionMatrix;
            
            model.gasContributionMatrix = model.compositionMatrix;
            model.gasContributionMatrix(:,~gasInd) = [];

            model.solidContributionMatrix = model.compositionMatrix;
            model.solidContributionMatrix(:,~solidInd) = [];
            
            model.compositionMatrix(:,phaseInd) = [];
            
            model.logSpeciesNames = cellfun(@(name) ['log', name], model.speciesNames, ...
                                         'uniformoutput', false);
                                     
            model.logGasNames = cellfun(@(name) ['log', name], model.gasNames, ...
                                         'uniformoutput', false);

            model.logSolidNames = cellfun(@(name) ['log', name], model.solidNames, ...
                                         'uniformoutput', false);
                                     
            model.nC = numel(model.speciesNames);
            model.nG = sum(gasInd);
            model.nS = sum(solidInd);
            
            for i = 1 : numel(model.elementNames)
                MCName = model.elementNames{i};
                ind = cellfun(@(x) ~isempty(x), regexpi(model.speciesNames,MCName));
                assert(sum(ind) > 0, ['The element ' MCName ' appears to make no contribution to a chemical species. Verify species and element names and/or remove the element if it is not needed.']); 
            end
            


        end

        %%
        function [model] = initChemicalReactions(model, reactions)
            %initChemicalReactions computes rection matrix as linear combinations of
            %species, must be called after initMasterComponents and
            %initSecondaryComponents
           
            assert(~isempty(model.speciesNames),'The function initChemicalReactions can only be called after initMasterComponents and initSecondaryComponents');

            % make sure input is a name value pair
            assert(rem(numel(reactions),2) == 0, 'Reaction and reaction constant must be given as a "key"/value pair.');

            rxn = regexprep(reactions(1:2:end),'[\s]','');

            rxnK = reactions(2:2:end);

            % unwrap master component vectors
            nRx = numel(rxn);
            nC = model.nC;
            nM = model.nMC;
            nG = model.nG;
            nS = model.nS;
            nP = model.nP;
            
            % check that an appropriate number of inputs have been
            % designated
            strlen = cellfun(@length,model.allComponentNames);
            [~, I] = sort(strlen, 2, 'descend');

            % make a copy of reactions
            rxnsO= rxn;

            % rewrite reactions in terms of the master component index
            % vectors
            for i = 1 : nRx;
                for j = 1 : nC + nS + nG + nP
                    ind = I(j);
                    [match,nomatch] = regexpi(rxn{i}, '(model\.Cind{\d})','match','split');
                    nomatch = strrep(lower(nomatch), lower(model.allComponentNames{ind}), ['model.Cind{' num2str(ind) '}']);
                    rxn{i} = strjoin(nomatch,match);
                end
                rxn{i} = strrep(rxn{i}, 'cind', 'Cind');
                remain = regexpi(rxn{i}, '[^(model\.Cind{\d})+-/*\.\d(<->)]','once');
                assert(isempty(remain),['The chemical reaction "' rxnsO{i} '" appears to contain a string combination that does not correspond to species. Ensure the species names contained in the reaction are consistent with those listed in secondaryComponents.']);
            end
            
            % check for <-> and handle it appropriately
            for i = 1 : nRx
                assert(~isempty(regexpi(rxn{i},'(=)','once')),...
                    ['Chemical reaction ' rxnsO{i} ' is improperly formatted. Reactions must be written as "reactants <-> products."']);
                rxn{i} = regexprep(rxn{i}, '=', ')+','ignorecase');
                try
                    tmp{i} = eval(['-(' rxn{i} ,';'])';
                catch
                    error(['The reaction "' rxnsO{i} '" appears to contain a string combination that does not correspond to species names. Make sure to add all relevant chemical species to the species list. If a reaction involves more than one species use "*" to designate this (i.e. 2*Na+ rather than 2Na+). ']);
                end
            end

            % assemble reaction matrix
            RM = vertcat(tmp{:});
            
            % assemble reaction constants
            model.reactionConstants = rxnK;
            model.logReactionConstants = cellfun(@(x) log(x), rxnK, 'UniformOutput', false);

            if any(cell2mat(cellfun(@(x) size(x), rxnK,'UniformOutput',false)) > 1)
                warning('A vector has been passed for the equilibrium constant of a reaction. Make sure the input vector size matches.')
            end
            % check that there are an equal number of each master component on both sides of the chemical reacitons.
            for i = 1 : nRx
                mCind = RM(i,:) ~= 0;
                mCcomp = repmat(RM(i,mCind),nM,1).*model.allContributionMatrix(:,mCind);
                mCcomp(strcmpi('e', model.elementNames),:) = 0;
                balance = sum(mCcomp,2);
                imbalanceInd = find(balance);
                for j = 1 : numel(imbalanceInd)
                    error(['The chemical reaction "' rxnsO{i} '" does not balance ' model.elementNames{imbalanceInd(j)} '.'])
                end
            end
            
            % check for charge balance within each reaction
            for i = 1 : nRx
                mCind = RM(i,:) ~= 0;
                mCcomp = RM(i,mCind).*model.allCharge(:,mCind);
                balance = sum(sum(mCcomp,2));
                assert(balance == 0, ['The chemical reaction "' rxnsO{i} '" does not balance charge.'])
            end
            
            model.nR = nRx;
            model.rxns = rxnsO;

            surfInfo = model.surfInfo;
                        
            SPMatrix = zeros(model.nR, nP);
            % Add surface activities to reaction matrix and remove
            % contribution to aqueous charge conservation
            if model.surfFlag && ~isempty(model.surfInfo)
                
                gNames = model.surfaces.groupNames;
                mNames = model.surfaces.masterNames;
                
                for i = 1 : numel(gNames)
                    for j = 1 : numel(mNames{i});
                        iterNames = mNames{i}{j};
                        mInd = strcmpi(iterNames, model.surfInfo.master);
                        sNames = surfInfo.species{mInd};
                        switch model.surfaces.scm{i}
                            case {'tlm','ccm'}
                                % find the number of layers
                                nL = numel(surfInfo.charge{mInd}{1});
                                layerInd = cellfun(@(x) ~isempty(x), regexpi(model.surfaceActivityCoefficientNames, [gNames{i} '_ePsi']));
                                for k = 1 : numel(sNames)
                                    % find the columns of the reaction matrix that contain
                                    % the species
                                    sp = sNames{k};
                                    spInd = strcmpi(sp, model.speciesNames);
                                    rmVec = RM(:,spInd);

                                    % replace the correct column and row with the
                                    % contribution
                                    SPMatrix(logical(rmVec),layerInd) = SPMatrix(logical(rmVec),layerInd) + repmat(rmVec(rmVec ~= 0),1,nL).*repmat(surfInfo.charge{mInd}{k},sum(logical(rmVec)),1);
                                end
                            case {'ie','langmuir'}
                                nL = 0;
                        end
                    end
                    numLay(i) = nL;
                end
            end

            % create reaction matrix

            model.allReactionMatrix = [RM SPMatrix];
    
            model.reactionMatrix = RM;
            model.reactionMatrix(:,model.phaseInd) = [];
            
            model.gasReactionMatrix = RM;
            model.gasReactionMatrix(:,~model.gasInd) = [];

            model.solidReactionMatrix = RM;
            model.solidReactionMatrix(:,~model.solidInd) = [];
            
            model.surfacePotentialMatrix = SPMatrix;
            
            assert(model.nR + model.nMC == model.nC + model.nG + model.nS, ['The given chemical system is rank deficient. The number of species must be equal to the sum of the number elements, surface functional groups, and reactions.'])
            
            nI = numel(model.chemicalInputModel.inputNames);
            
            assert(nI == model.nMC , ['For the defined chemical system ' num2str(model.nMC) ' species, elements or species, must be designated as inputs. Use an "*" to designated a component as an input.']);


            model.allContributionMatrix = [model.allContributionMatrix zeros(model.nMC, nP)];
            
            % make sure multidentate species do not span multiple surfaces
            if model.surfFlag && ~isempty(model.surfInfo)
                first = 1;
                for i = 1 : numel(gNames)
                    nL = numLay(i);
                    psiPack = SPMatrix(:,first:first+nL-1);
                    first = first + nL;
                    stack(:,i) = logical(sum(psiPack,2));
                end

                mdInd = sum(stack,2) > 1;
%                 if any(mdInd)
%                     error('Multidentate species may not span multiple surface groups.') 
%                 end
            end
        end

        %%
        function model = setupMaxMatrices(model)

            C = model.compositionMatrix;
            nC = model.nC;
            maxMatrices = cell(nC, 1);
            for i = 1 : nC
                maxMatrices{i} = diag(1./C(:, i));
            end
            model.maxMatrices = maxMatrices;

        end


        %%
        function model = initLinearCombinations(model, p)
            
            linearCombos = p.Results.combinations;
            
            if isempty(linearCombos)
               	model.nLC = 0;
                model.combinationNames = {};
                return
            end
            
            %initElectrostaticModel parses the surface chemistry inputs 
            
            newVar = linearCombos(1:2:end);
            combos = linearCombos(2:2:end);
            
            nVar = numel(newVar);
            
            % find input components and remove the asterik from the name
            InInd = cellfun(@(x) ~isempty(x), regexp(newVar, '*'));
            names = regexprep(newVar, '*','');



            if isempty(model.chemicalInputModel)
                model.chemicalInputModel = ChemicalInputModel();
            end
            model.chemicalInputModel.inputNames = ...
                horzcat(model.chemicalInputModel.inputNames, names(InInd));
            
            
            % unwrap master component vectors
            nS = numel(model.allComponentNames);
            nM = model.nMC;

            % sort comp names by length
            strlen = cellfun(@length,model.allComponentNames);
            [~, I] = sort(strlen, 2, 'descend');
            
            combos0 = combos;
            
            phaseNames = [model.gasNames, model.solidNames];
            for i = 1 : numel(phaseNames)
                testFlag = any(cellfun(@(x) ~isempty(x), strfind(combos, phaseNames{i}), 'uniformoutput', true));
                if testFlag
                    error('Pure phases can not be used in linear combinations.');
                end
            end
            
            % rewrite combinations in terms of the master component index
            % vectors
            for i = 1 : nVar;
                for j = 1 : nS
                    ind = I(j);
                    [match,nomatch] = regexpi(combos{i}, '(model\.Cind{\d})','match','split');
                    nomatch = strrep(lower(nomatch), lower(model.allComponentNames{ind}), ['model.Cind{' num2str(ind) '}']);
                    combos{i} = strjoin(nomatch,match);
                end
                combos{i} = strrep(combos{i}, 'cind', 'Cind');
            end 
            
            for i = 1 : nVar
                 
                try
                    tmp{i} = eval([ combos{i} ,';'])';
                    if numel(tmp{i}) == 1
                        tmp{i} = zeros(1, model.nC);
                    end
                catch
                    error(['The linear combination "' combos0{i} '" appears ' ...
                                        'to contain a string combination that ' ...
                                        'does not correspond to species names. ' ...
                                        'Make sure to add all relevant chemical ' ...
                                        'species to the species list. Use  ' ...
                                        '"*" after the coefficients in the ' ...
                                        'linear combinations  (i.e. 2*Na+ rather than 2Na+). ']);
                end
            end
            
            
            % create combination matrix
            model.allCombinationMatrix = vertcat(tmp{:});
            model.combinationMatrix = model.allCombinationMatrix;
            model.combinationMatrix(:, model.phaseInd) =[];
            
            model.combinationNames = names;
            model.nLC = numel(names);
            
           

        end
%%
        function model = initElectrostaticModel(model, givenSurfaceInformation, varargin)
            %initElectrostaticModel parses the surface chemistry inputs 
            
            if isempty(givenSurfaceInformation)
                model.nP = 0;
                model.surfaceActivityCoefficientNames = cell(1,0);
                model.logSurfaceActivityCoefficientNames = cell(1,0);
                return
            end
            model.surfaceActivityCoefficientNames = {};
            
            funcGroup = givenSurfaceInformation(1:2:end);
            addInfo = givenSurfaceInformation(2:2:end);
            
            funcGroup = regexprep(funcGroup,'[\s]','','ignorecase');

            % some checks to get the input format correct
            assert(iscell(addInfo), 'Surface information must be a cell.');
                        
            surfInfo.master = funcGroup;

            % grab surface area and slurry density inputs
            try
                surfInfo.d = cellfun(@(x) x{1}(:,1), addInfo, 'UniformOutput', false);
                surfInfo.s = cellfun(@(x) x{1}(:,2), addInfo, 'UniformOutput', false);
                surfInfo.a = cellfun(@(x) x{1}(:,3), addInfo, 'UniformOutput', false);
            catch
                % better error checking here
                error('Surface site density, specific surface area, and slurry density must be given as an array of doubles in that order.')
            end

            % identify electrostatic model
            options = {'tlm', 'ccm', 'langmuir','ie'};
            modNames = cellfun(@(x) x{2}, addInfo,'UniformOutput', false);
            surfInfo.scm = cellfun(@(x) validatestring(x,options,'initElectrostaticModel','electrostatic model'), modNames, 'UniformOutput',false);

            % make sure correct input size is provided for the surface model
            chargeString = 'The charge of each surface species for this type of model is taken from the species names.';
            for i = 1 : numel(surfInfo.master) 
               	t = numel(addInfo{i}(3:end));
                
                switch surfInfo.scm{i}
                	case {'langmuir', 'ie'}
                        assert(t == 0, ['Too many input arguments have been given to the ' surfInfo.master{i} ' informtion cell. For the ' surfInfo.scm{i} ' model only the surface geometry and the designation "' surfInfo.scm{i} '" are needed. ' chargeString]);
                    case 'ccm'
                        if t == 0
                            error(['Too few arguments have been given to the ' surfInfo.master{i} ' informtion cell. As a constant capacitance surface the geometry, "ccm" tag, and capacitance value must be provided.']);
                        elseif t > 1
                            error(['Too many arguments have been given to the ' surfInfo.master{i} ' informtion cell. As a constant capacitance surface only the geometry, "ccm" tag, and capacitance value must be provided. ' chargeString]);   
                        end
                    case 'tlm'
                      	assert(t >= 1, ['Too few arguments have been given to the ' surfInfo.master{i} ' informtion cell. As a triple layer surface the geometry, "tlm" designation, and capacitance values must be provided.']);
                           
                end
            end                
            
            % grab capacitance vales
            for i = 1 : numel(surfInfo.scm)
                switch surfInfo.scm{i}
                    case 'ccm'
                        surfInfo.c{i} = addInfo{i}{3};
                        assert(size(surfInfo.c{i},2) == 1, 'If the constant capacitance model is specified only one value for capacitance is expected.')
                    case 'tlm' 
                        surfInfo.c{i} = addInfo{i}{3};
                        assert(size(surfInfo.c{i},2) == 2, 'If the triple layer model is specified 2 values for capacitance are expected. To simulate the diffuse layer model provide a high capacitance to both entries. To simulate the basic stern model provide a high capacitance value to the second layer.');
                end
            end
    
            for i = 1 : numel(surfInfo.master)
                MName = surfInfo.master{i};
                
                mInd = strcmpi(MName, model.elementNames);
                assert(any(mInd), ['Given surface master species ', MName ' does not match known master species']);
                ind = logical(model.compositionMatrix(mInd,:));
                speciesNames{i} = model.speciesNames(ind);
                
                
                % grab the charge from the species name
                nS = numel(speciesNames{i});
                tmpvec = zeros(1, nS);
                
                for j = 1 : nS
                    sInd = strcmpi(speciesNames{i}{j}, model.speciesNames);
                    tmp = model.allCharge(sInd);
                    charge{i}{j} = tmp;
                    namedCharge{i}{j} = tmp;
                end
                
                givenNames = addInfo{i}(4:2:end);
                givenNames = regexprep(givenNames,'[\s]','','ignorecase');

                givenCharge = addInfo{i}(5:2:end);
                
                assert(numel(givenNames) == numel(givenCharge), 'If a name of a species is given to the surfaces cell, a charge value is expected.');
                
                % make sure ion exchange surfaces make sense
                if strcmpi(surfInfo.scm{i}, 'ie')
                    assert(all(cellfun(@(x) x == 0, charge{i})), 'A charge has been assigned to an ion exchange surface species. This is contradictory to the chemical model. Consider using a different model, or remove the charge designation.')
                end
                
                switch surfInfo.scm{i}
                    case 'tlm'
                        for j = 1 : nS
                            spInd = strcmpi(speciesNames{i}{j},givenNames);
                            if any(spInd)
                                charge{i}{j} = givenCharge{spInd}(:)';
                            end
                            charge{i}{j} = horzcat( charge{i}{j}, zeros(1, 3 - numel(charge{i}{j})));
                            assert(sum(charge{i}{j}) == namedCharge{i}{j}, ['The charge distribution of the surface species ' speciesNames{i}{j} ' does not sum to the total species charge as determined by the species name.']);

                        end  

                        for j = 1 : numel(givenNames)
                            assert(logical(sum(strcmpi(givenNames{j}, speciesNames{i}))), ['The surface species ' givenNames{j} ' is not tabulated in species.']);
                        end
                end
                
            end
            surfInfo.species = speciesNames;
            surfInfo.charge = charge;
            
            model.surfInfo = surfInfo;

            % add surface activities to composition matrix
            for i = 1 : numel(model.surfaces.groupNames)
                model.surfaceChargeNames = {};
                model.surfacePotentialNames = {};
                m = model.surfaces.groupNames{i};
                
                mastNames = model.surfaces.masterNames{i};

                
                % make sure model type and capacitances are the same for
                % all surfaces in a group
                for j = 1 : numel(mastNames)
                    testInd = strcmpi(mastNames{j},model.surfInfo.master);
                    mTest{j} = model.surfInfo.scm{testInd};
                    cTest{j} = 0;
                    switch mTest{j}
                        case {'ccm','tlm'}
                            cTest{j} = model.surfInfo.c{testInd};
                    end
                    
                end
                
                assert(isequal(mTest{:},mTest{:}), 'Surface functional groups that are combined into a single surface must use the same electrostatic model.');
                assert(isequal(cTest{:},cTest{:}), 'Surface functional groups that are combined into a single surface must use the same values for capacitance.');

                model.surfaces.scm{i} = mTest{1};
                model.surfaces.c{i} = cTest{1};                

                switch mTest{1}
                    case 'ccm'
                        model.allComponentNames{end+1} = [ m '_ePsi'];
                        model.surfaceActivityCoefficientNames{end+1} = [ m '_ePsi'];
                        model.surfaceChargeNames{end+1} = [m '_sig'];
                        model.surfacePotentialNames{end+1} = [m '_Psi'];
                    case 'tlm'
                        model.allComponentNames =  horzcat( model.allComponentNames, [m '_ePsi_0'], [m '_ePsi_1'], [m '_ePsi_2']);
                        model.surfaceActivityCoefficientNames =  horzcat( model.surfaceActivityCoefficientNames, [m '_ePsi_0'], [m '_ePsi_1'], [m '_ePsi_2']);
                        model.surfaceChargeNames = horzcat(model.surfaceChargeNames, [m '_sig_0'],[m '_sig_1'],[m '_sig_2']);
                        model.surfacePotentialNames = horzcat(model.surfacePotentialNames, [m '_Psi_0'],[m '_Psi_1'],[m '_Psi_2']);
                end
            end

            model.logSurfaceActivityCoefficientNames = cellfun(@(name) ['log', name], model.surfaceActivityCoefficientNames, ...
                                         'uniformoutput', false);

            model.nP = numel(model.surfaceActivityCoefficientNames);


        end
        
        
        %%
        function printChemicalSystem(model, varargin)
            % printChemicalSystem prints the mass conservation, law of mass
            % action, and linear combination equations to the workspace 
            
            chemicalSystemPrintFunction(model, varargin{:});
          
        end


        %%
        function [state, report] = updateState(model, state, problem, dx, drivingForces) %#ok
        % Update state based on Newton increments
            state0 = state;
            [state, report] = updateState@PhysicalModel(model, state, problem, ...
                                                        dx, drivingForces);

            state = model.syncFromLog(state);
            state = updateChemicalModel(model, problem, state, state0 );
            state = model.syncLog(state);
            

            if model.plotIterations
                h = findobj('tag', 'updatechemfig');
                if isempty(h)
                    figure
                    set(gcf, 'tag', 'updatechemfig');
                    h = findobj('tag', 'updatechemfig');
                end
                set(0, 'currentfigure', h)
                clf
                if size(state.species, 1) ~= 1
                    plot(log10(state.species*litre/mol));
                    title('species (chemistry step)');
                    legend(model.speciesNames);
                end

                h = findobj('tag', 'chemmastercomp');
                if isempty(h)
                    figure
                    set(gcf, 'tag', 'chemmastercomp');
                    h = findobj('tag', 'chemmastercomp');
                end
                set(0, 'currentfigure', h)
                clf
                if size(state.species, 1) ~= 1
                    plot(log10(state.elements*litre/mol));
                    title('master species (chemistry step)');
                    legend(model.elementNames);
                end
                drawnow;
            end

        end

        %%
        function state = validateState(model, state)
            state = validateState@PhysicalModel(model, state);
            state = model.syncLog(state);
        end

        %%
        function state = syncLog(model, state)
            state.logElements = log(state.elements);
            state.logSpecies       = log(state.species);
            
            if isfield(state, 'saturationIndicies');
                state.logSaturationIndicies     = log(state.saturationIndicies);
            end
            if isfield(state, 'partialPressures');
                state.logPartialPressures       = log(state.partialPressures);
            end
            if isfield(state, 'surfaceActivityCoefficients');
                state.logSurfaceActivityCoefficients       = log(state.surfaceActivityCoefficients);
            end
        end

        %%
        function state = syncFromLog(model, state)
            state.elements = exp(state.logElements);
            state.species       = exp(state.logSpecies);
            
            if isfield(state, 'logSaturationIndicies');
                state.saturationIndicies     = exp(state.logSaturationIndicies);
            end
            if isfield(state, 'logPartialPressures');
                state.partialPressures       = exp(state.logPartialPressures);
            end
            if isfield(state, 'logSurfaceActivityCoefficients');
                state.surfaceActivityCoefficients       = exp(state.logSurfaceActivityCoefficients);
            end
        end


        %%
        function [logComps, logMasterComps] = prepStateForEquations(model, state)

            logComps = cell(model.nC, 1);
            [logComps{:}] = model.getProps(state, model.logSpeciesNames{:});

            logMasterComps = cell(model.nMC, 1);
            [logMasterComps{:}] = model.getProps(state, model.logElementNames{:});

            [logComps{:}] = initVariablesADI(logComps{:});

        end

        %%
        function [state, failure, report] = solveChemicalState(model, inputstate)
        % inputstate contains the input and the initial guess.

            inputstate = model.validateState(inputstate); % in particular,
                                                          % updates the log
                                                          % variables if necessary.

            solver = NonLinearSolver('maxIterations', model.maxIterations);
            dt = 0; % dummy timestep
            drivingForces = []; % drivingForces;
            inputstate0 = inputstate;

            [state, failure, report] = solveMinistep(solver, model, inputstate, ...
                                                     inputstate0, dt, ...
                                                     drivingForces);

        end

        %%
        function [state, model] = computeActivities(model, state)
        % computeAcitivities computes the acitivity of each aqueous species
        % using the extended Davies equaiton and adds the values to the 
        % field state.activities.
        %
        % SYNOPSIS:
        %  [state, chem] = computeActivities(chem, state)
        %
        %
        % REQUIRED PARAMETERS:
        %   state - the state variable produced by chem.initState. 
        %          
        %
        % RETURNS:
        %   state - The state structure containing the additional field 
        %           state.acitivities which holds the acitivity of each aqueous species 
        %           in units of mol/meter^3. Activities can be retrieved using the 
        %           getProp command by calling for the species name prepended by 'a.'
        %
        %   chem - updated chemical object including the field chem.activityNames 
        %          for use in calling the activities using getProp/s.
        %
        % EXAMPLE:
        %   [state, chem] = chem.computeActivities(state);
        %   aH2O = chem.getProp(state, 'aH2O');
        %   pH = -log10(chem.getProp, 'aH+');
        %
        % SEE ALSO:
        %   'computeChargeBalance'
        
        %{
        Copyright 2009-2017 SINTEF Digital, Mathematics & Cybernetics.

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
        
            assert( nargout == 2, ['Output argument to computeActivities must include the state variable and the chemical model object']);
            [state, model] = activity(model, state);
            
        end
        
        %%
        function [state, model] = computeChargeBalance(model, state)
        % computeChargeBalance computes the residual of the aqueous charge
        % balance equation and adds the values to the field state.chargeBalance.
        %
        % SYNOPSIS:
        %  [state, chem] = computeChargeBalance(chem, state)
        %
        %
        % REQUIRED PARAMETERS:
        %   state - the state variable produced by chem.initState.
        %
        % RETURNS:
        %   state - The state structure containing additional field state.chargeBalance
        %           which holds the aqueous charge balance residual in units of fraction of total
        %           charge. The charge balance can be retrieved using the getProps
        %           command by calling for the variables 'chargeBalance'.
        %
        %   chem  - updated chemical object including for use in calling the
        %           charge balance residual using getProp/s.
        %
        % EXAMPLE:
        %   [state, chem] = chem.computeChargeBalance(state);
        %   charge = chem.getProp(state, 'chargeBalance');
        %
        % SEE ALSO:
        %   'computeActivities'
        
        %{
        Copyright 2009-2017 SINTEF Digital, Mathematics & Cybernetics.

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

            assert( nargout == 2, ['Output argument to computeChargeBalance must include the state variable and the chemical model object']);
            [state, model] = chargeBalance(model, state);
            
        end
        
        %%
        function [state, model] = computeSurfacePotentials(model, state)
        % computeSurfacePotentials computes the potential of each layer of each
        % surface and adds the values to the field state.surfacePotentials.
        %
        % SYNOPSIS:
        %  [state, chem] = computeSurfacePotentials(chem, state)
        %
        % REQUIRED PARAMETERS:
        %   state - the state variable produced by model.initState.
        %          
        % RETURNS:
        %   state - The state structure containing the additional field state.surfacePotentials of
        %           which holds the potential of each layer of
        %           each surface in Volts. Values can be retrieved using
        %           the getProps command by calling for the surface functional
        %           group name, followed by '_Psi_' followed by the layer
        %           number (0, 1, 2) or for the constant capacitance model just
        %           use '_Psi'.
        %
        %   chem - updated chemical object including the field chem.surfacePotentialNames 
        %           for use in calling the surface potentials using getProp/s.
        %
        % EXAMPLE:
        %   [state, chem] = chem.computeSurfacePotentials(state);
        %   [pot1, pot2] = chem.getProps(state, '>SiO_Psi_0', '>SiO_Psi_1');
        %
        % SEE ALSO:
        %   'computeSurfaceCharges'
        
        %{
        Copyright 2009-2017 SINTEF Digital, Mathematics & Cybernetics.

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
        
            assert( nargout == 2, ['Output argument to computeSurfacePotentials must include the state variable and the chemical model object']);
            [state, model] = surfacePotential(model, state);
            
        end
        
        %%
        function [state, model] = computeSurfaceCharges(model, state)
        % computeSurfaceCharges computes the charge of each layer of each
        % surface and adds the values to the field state.surfaceCharges.
        %
        % SYNOPSIS:
        %  [state, chem] = computeSurfaceCharges(chem, state)
        %
        % REQUIRED PARAMETERS:
        %   state - the state variable produced by chem.initState.
        %          
        % RETURNS:
        %   state - The state structure with the additional field 
        %           state.surfaceCharges containing the charge of each layer of
        %           each surface in C/meter^2. Values can be retrieved using
        %           the getProps command by calling for the surface functional
        %           group name, followed by '_sig_' followed by the layer
        %           number (0, 1, 2) or for the constant capacitance model just
        %           use '_sig'.
        %
        %   chem - updated chemical object including the field chem.surfaceChargeNames 
        %           for use in calling the surface charges using getProp/s.
        %
        % EXAMPLE:
        %   [state, chem] = chem.computeSurfaceCharges(state);
        %   charge0 = chem.getProp(state, '>SiO_sig_0');
        %
        % SEE ALSO:
        %   'computeSurfacePotentials'
        
        %{
        Copyright 2009-2017 SINTEF Digital, Mathematics & Cybernetics.

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
            assert( nargout == 2, ['Output argument to computeSurfaceCharges must include the state variable and the chemical model object']);                
            [state, model] = surfaceCharge(model, state);
        end
        
        %%
        function [state, model] = computeAqueousConcentrations(model, state)
        % computeAqueousConcentrations computes the total concentration of each 
        % element in the aqueous phase. This excludes, gas, solid and
        % surface bound species. 
        %
        % SYNOPSIS:
        %  [state, chem] = computeAqueousConcentrations(chem, state)
        %
        % REQUIRED PARAMETERS:
        %   state - the state variable produced by chem.initState.
        %          
        % OUTPUTS:
        %   state - The state structure containing the additional field 
        %           state.aqueousConcentrations which hold the total aqueous
        %           concentration of each element in the aqueous phase. The values can be retrieved using the
        %           getProps function by asking for the element name followed
        %           by '(aq)'.
        %
        %   chem - updated chemical object including the field chem.surfaceChargeNames 
        %          for use in calling the surface charges using getProp/s.
        %
        % EXAMPLE:
        %   [state, chem] = chem.computeAqueousConcentrations(state);
        %   Na = chem.getProp(state, 'Na(aq)');
        %   H = chem.getProp(state, 'H(aq)');
        %
        % SEE ALSO:
        %   'computeSurfaceConcentrations'
        
        %{
        Copyright 2009-2017 SINTEF Digital, Mathematics & Cybernetics.

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
        
            assert( nargout == 2, ['Output argument to computeAqueousConcentrations must include the state variable and the chemical model object']);
            [state, model] = aqueousConcentrations(model, state);
            
        end
        
        %%
        function [state, model] = computeSurfaceConcentrations(model, state)
        % computeSurfaceConcentrations computes the total concentration of each 
        % element bound to surfaces.
        %
        % SYNOPSIS:
        %  [state, chem] = computeAqueousConcentrations(chem, state)
        %
        % REQUIRED PARAMETERS:
        %   state - the state variable produced by model.initState.
        %          
        % OUTPUTS:
        %   state - The state structure containing the additional field 
        %           state.surfaceConcentrations which holds the surface
        %           concentration of each element. The values can be retrieved using the
        %           getProps function by asking for the element name followed
        %           by '_surf'.
        %
        %   chem - updated chemical object now includes the field 
        %           chem.surfaceConcentractionNames for use in  calling the 
        %           surface concetrations using getProp/s. 
        %
        % EXAMPLE:
        %   [state, chem] = chem.computeSurfaceConcentrations(state);
        %   [Na, H] = chem.getProps(state, 'Na(surf)', 'H(surf)');
        %
        % SEE ALSO:
        %   'computeSurfaceConcentrations'
        
        %{
        Copyright 2009-2017 SINTEF Digital, Mathematics & Cybernetics.

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
        
            assert( nargout == 2, ['Output argument to computeSurfaceConcentrations must include the state variable and the chemical model object']);
            [state, model] = surfaceConcentrations(model, state);
            
        end
    end

end
