classdef ChemicalModel < PhysicalModel
 %ChemicalModel A model for solving equilibrium geochemistry.
 
 %{
            Copyright 2009-2017 SINTEF DIGITAL, Applied Mathematics and Cybernetics and The University of Texas at Austin.

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
        MasterCompNames      % Names of the master components as first master component.
        logMasterCompNames   % Names of the log master components (= 'log' in
                             % front of MasterCompNames).

        MCind                % index vector of each master component
        nMC                  % number of master components

        CompNames            % Names of the components.
        logCompNames         % Names of the log components (= 'log' in
                             % front of CompNames).

        CompositionMatrix    % Composition for each components, in term of the
                             % master components.
        Cind                 % index vector of secondary components
        nC                   % Number of components
        
        maxMatrices          % cells of matrices, one per components. Used to
                             % compute the upper bound when updating a component.

        ChargeVector         % coefficients of charged species.
        ReactionMatrix       % Coefficients for the reaction.
        ReactionConstants    % Contains the K_i.
        LogReactionConstants % Contains the log K_i.


        nR                   % Number of reactions
        rxns                 % string of reactions

        maxIterations        % max nonlinear iterations for solveChemicalState

        surfChargeMatrix     % matrix of surface charge contributions
        surfMaster           % index vector of master components that are surfaces

        surfFlag             % are there surface species 1 or 0
        surfaces             % surface groups for electrostatic calculations
        
        nSurf                % number of surface components
        
        compositionModel     % model for computing the mass conservation solution
        compositionReactionModel % model for computing the mas conservation and reaction equations together
        chargeBalanceModel   % model for computing the solution when charge balance is required
        
        chemicalInputModel   % Model to solve chemical equilibrium state from user
                              % input. see initState member function. variable
                              % is initialized in initSecondaryComponents.

        surfInfo            % information regarding the surfaces
        CompActivityNames   % name of activity of species
        
        plotIter            % toggle plotting of iteration information true or false
        
        surfaceChargeNames         % name of surface charge variables
        surfacePotentialNames         % name of potentials
        
        allCharge           % vector for ensuring each reaciton is charge balanced
        
        inputs              % names of inputs the user must provide
        
        CombinationNames    % names of linear combinations
        CombinationMatrix   % combination for each linear combination in terms of secondary components
        nLC                 % number of linear combinations
        
        AqueousConcentrationNames % name of total aqueous concentration variables
        SurfaceConcentrationNames % name of total surface concentration variables
        
        SolidNames          % name of species in the solid phase
        GasNames            % name of species in the gas phase
        
    end


    methods

        %%
        function model = ChemicalModel(varargin)
        %ChemicalModel A model for solving equilibrium geochemistry.
            %
            % SYNOPSIS:
            %  chem = ChemicalModel(varargin)
            %
            % DESCRIPTION:
            %   A class of PhysicalModel which can construct and solve
            %   arbitrary aqueous geochemical system including surface chemistry
            %   all under the assumption of local chemical equilibrium. 
            %
            % REQUIRED PARAMETERS:
            %   elementNames        - A cell array of strings containing all
            %       elements to be considered in the chemical system. 
            %       Elements do not have to correspond to actual element names,
            %       but it is reccomened that they do.
            %            
            %           elementNames = {'H', 'O'};
            %
            %   speciesNames        - A cell array of strings of the chemical
            %       species to be considered in the chemical system. Each 
            %       entry in speciesNames must be a combination of
            %       entries found in elementNames, or surfaces. The species
            %       charge can be desginated with a +/- at the end of the
            %       name followed by the charge number (i.e. Ca+2, Cl-).
            %
            %           speciesNames = {'H+', 'OH-', 'H2O'};
            %
            %   reactions           - A cell array listing the chemical
            %       reactions of the chemical system as strings, followed
            %       by the equilibrium constant of the reaction as a scalar
            %       quantity. Entries in reactions must be given as a
            %       string in the form 'reactants <-> products'. 
            %       Equuilibrium constant must be given in SI units.
            %
            %           reactions = {'H2O <-> H+ + OH-', 1e-14*mol/litre};
            %
            % OPTIONAL PARAMETERS:
            %   surfaces            - A cell structure containing the names
            %       of surface functional groups, followed by information
            %       regarding the specific surface. The first entry of
            %       surfaces is the name of the surface functional group.
            %       Surface functional group names must begin with the ">"
            %       symbol (i.e. '>FeO'). These names shoudl be used in the same
            %       way an entries in elementNames for the construction of
            %       speciesNames. There are three different surface
            %       chemistry models implemented in this ChemicalModel: the
            %       langmuir model, the triple layer model, and the
            %       constant capacitance model. The diffuse layer and basic
            %       stern models can be constructed by an adjustment of the
            %       parameters of the triple layer model. A chemical system
            %       can have any number of surfaces with different surface
            %       complexation models.
            %
            %       The Langmuir model is the most simplistic case.
            %
            %           geometry = [1/(nano*meter)^2 50*meter^2/gram 1000*grams/litre];
            %
            %           surfaces = {'>FeO', {geometry, 'langmuir'}
            %
            %       In this example the surface functional group is
            %    	'>FeO'. The second entry is a cell containing the
            %     	surface parameters. The variable geometry must
            %   	always be the first entry in this parameters cell.
            %     	The entries of geometry must correspond to the
            %      	surface site density, the specific surface area, and
            %      	the slurry density in that specific order. Values
            %     	must be given in SI units. The second entry of the
            %     	parameters cell must be the type of surface model to
            %     	use for the surface. In the case of 'langmuir' as above
            %     	no other entries are needed. 
            %
            %       The triple layer model requires additional work.
            %
            %           capacitances = [1 0.2]*Farad/meter^2;
            %
            %           surfaces = {'>FeO', {geometry, 'tlm', capacitance,...
            %                                           '>FeO-', [-1 0 0]}
            %
            %       In this example geometry is defined as before, but now
            %       the '>FeO' is a triple layer model surface, as is
            %       indicated by 'tlm.' Additional parameters are need in
            %       this case. The capacitance density of each of the
            %       Helmholtz layers must be given after the model type.
            %       The basic stern model can be simulated by increasing
            %       the capacitance of the second layer (i.e. [1 1e3]) and
            %       the diffuse layer model can be simulated by increasing
            %       the capacitance of both layers (i.e. [1e3 1e3]). The
            %       capacitance is followed by the name of each surface
            %       species associated with '>FeO' and the species charge
            %       contribution to each plane. 
            %
            %       The constant capacitance model can be similarly
            %       defined.
            %
            %           capacitances = 1*Farad/meter^2;
            %
            %           surfaces = {'>FeO', {geometry, 'ccm', capacitance, '>FeO-', -1}
            %
            %       In this example geometry is defined as before, but now
            %       the '>FeO' is a constant capacitance surface, as is
            %       indicated by 'ccm.' Additional parameters are need in
            %       this case. The capacitance density of the
            %       Helmholtz layers must be given after the model type.
            %       The capacitance is followed by the name of each surface
            %       species associated with '>FeO' and the species charge
            %       contribution to the surface.       
            
            %{
            Copyright 2009-2017 SINTEF DIGITAL, Applied Mathematics and Cybernetics and The University of Texas at Austin.

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
            
            model = model@PhysicalModel([]);
%             model.maxIterations = 50;
            if nargin > 3
                % Creating ChemicalModel
                MasterCompNames = varargin{1};
                CompNames       = varargin{2};
                Reactions       = varargin{3};
                
            	p = inputParser;
            
                valFun = @(x) iscell(x);
                p.addOptional('surfaces', '', valFun);                
                p.addOptional('combinations', '', valFun);
                p.addOptional('electrostaticgroup', '', valFun);

                p.parse(varargin{4:end})
            
                % add surface functional groups to master component list
                if ~isempty(p.Results.surfaces)
                    
                    assert(rem(numel(p.Results.surfaces),2)==0, 'A key/value pair must be provided for surfaces, yet the number of inputs is odd.');
                    
                    surfNames = p.Results.surfaces(1:2:end);
                    masterNames = surfNames;
                    
                    toPass = p.Results.surfaces;
                    
%                     ind = cellfun(@(x) ~isempty(x), regexpi(surfNames,'g[roup]'));
%                     masterNames = surfNames(~ind);
                    
                    if ~isempty(p.Results.electrostaticgroup)
%                         tmp = find(ind);

                        % seperate group names and master components that
                        % comprise the groups
                        groups = p.Results.electrostaticgroup;
                        groupNames = groups(:,1);
                        groupMasterNames = groups(:,2:end);
                        
%                         toPass = p.Results.surfaces;
%                         toPass(2*tmp) =[];
%                         toPass(2*tmp-1) = [];
                        
                        % make sure group names are unique
                        for i = 1 : numel(groupNames)
                            ind = zeros(size(groupNames));
                            ind(i) = 1;
                            test = true;
                            if numel(groupNames) == 2
                                test = sum(strcmpi(groupNames{i}, groupNames{i+1}));
                                test = test == 0;
                            elseif numel(groupNames) > 2
                                test = all(isempty(strcmpi(groupNames{i}, groupNames{~ind})));
                            end
                            assert(test, 'Surface group names must be unique.');
                        end
                        
                        % make sure surface master species do not appear in
                        % more than one group
                        for i = 1 : numel(groupMasterNames)
                            ind = zeros(size(groupMasterNames));
                            ind(i) = 1;
                            test = true;
                            
                            if numel(groupNames) == 2
                                test = sum(strcmpi(groupMasterNames{i}, groupMasterNames{i+1}));
                                test = test == 0;
                            elseif numel(groupMasterNames) > 2
                             	test = ~all(strcmpi(groupMasterNames{i}, groupMasterNames(~ind)));
                            end
                            
                            validatestring(groupMasterNames{i}, surfNames, 'ChemicalModel', 'entries to share');
                            assert(test, ['The surface species ' groupMasterNames{i} ' is repeated or appears in more than one group.']);
                        end
                        
                        % add functional groups not listed in groups to
                        % groupNames may need to change cell indexing for when there
                        % arent the same number of surface master components in  a group
                        for i = 1 : numel(surfNames);
                            if ~any(strcmpi(surfNames{i}, groupMasterNames))
                                groupNames{end+1} = surfNames{i};
                                groupMasterNames{end+1,1} = surfNames{i};
                            end
                        end
                    else
                        
                        groupNames = surfNames';
                        groupMasterNames = surfNames';
                    end
                    
                    model.surfaces.groupNames = groupNames;
                    model.surfaces.groupMasterNames = groupMasterNames;
                    
                    surfNames = cellfun(@(name) [name '*'], ...
                                               surfNames, ...
                                               'uniformoutput', false);
                                           
                    MasterCompNames = horzcat(MasterCompNames, surfNames);
                end
                
                model = initMasterComponents(model, MasterCompNames);
                model = initSecondaryComponents(model,  CompNames);
                
                if ~isempty(p.Results.surfaces)
                    model = initElectrostaticModel(model, toPass);
                end
                
                if ~isempty(p.Results.combinations)
                    model = initLinearCombinations(model, p.Results.combinations);
                else
                    model.nLC = 0;
                    model.CombinationNames = {};
                end
                
                model = initChemicalReactions(model, Reactions{:});
                model = setupMaxMatrices(model);
                props = properties(model);
                for i = 1 : numel(props);
                    if (~strcmpi('compositionModel', props{i})) && (~strcmpi('compositionReactionModel', props{i}))
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

                model.plotIter = false;
        end

        %%
        function [problem, state] = getEquations(model, state0, state, dt, drivingForces, varargin)

            [logcomps, logmasterComps] = prepStateForEquations(model, state);

            [eqs, names, types] = equationsChemicalLog(logcomps, logmasterComps, ...
                                                       model);

            primaryVariables = model.logCompNames;
            problem = LinearizedProblem(eqs, types, names, primaryVariables, state, dt);

        end



        %%
        function [fn, index] = getVariableField(model, name)

            varfound = false;

            while ~varfound

                if strcmpi(name, 'components')
                    varfound = true;
                    fn = 'components';
                    index = ':';
                    break
                end

                if strcmpi(name, 'chargeBalance')
                    varfound = true;
                    fn = 'chargeBalance';
                    index = ':';
                    break
                end
                
                if strcmpi(name, 'masterComponents')
                    varfound = true;
                    fn = 'masterComponents';
                    index = ':';
                    break
                end

                if strcmpi(name, 'logcomponents')
                    varfound = true;
                    fn = 'logcomponents';
                    index = ':';
                    break
                end
                
                if strcmpi(name, 'activities')
                    varfound = true;
                    fn = 'activities';
                    index = ':';
                    break
                end

                if strcmpi(name, 'logmasterComponents')
                    varfound = true;
                    fn = 'logmasterComponents';
                    index = ':';
                    break
                end

                if strcmpi(name, 'combinationComponents')
                    varfound = true;
                    fn = 'combinationComponents';
                    index = ':';
                    break
                end
                 
                ind = strcmpi(name, model.CombinationNames);
                if any(ind)
                    varfound = true;
                    fn = 'combinationComponents';
                    index = find(ind);
                    break
                end
              
                
                ind = strcmpi(name, model.CompNames);
                if any(ind)
                    varfound = true;
                    fn = 'components';
                    index = find(ind);
                    break
                end
                
                if strcmpi(name, 'aqueousConcentrations')
                    varfound = true;
                    fn = 'aqueousConcentrations';
                    index = ':';
                    break
                end

                ind = strcmpi(name, model.AqueousConcentrationNames);
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

                ind = strcmpi(name, model.SurfaceConcentrationNames);
                if any(ind)
                    varfound = true;
                    fn = 'surfaceConcentrations';
                    index = find(ind);
                    break
                end
                
                ind = strcmpi(name, model.CompActivityNames);
                if any(ind)
                    varfound = true;
                    fn = 'activities';
                    index = find(ind);
                    break
                end
                
                ind = strcmpi(name, model.logCompNames);
                if any(ind)
                    varfound = true;
                    fn = 'logcomponents';
                    index = find(ind);
                    break
                end
                
                ind = strcmpi(name, model.MasterCompNames);
                if any(ind)
                    varfound = true;
                    fn = 'masterComponents';
                    index = find(ind);
                    break
                end

                ind = strcmpi(name, model.logMasterCompNames);
                if any(ind)
                    varfound = true;
                    fn = 'logmasterComponents';
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
                
                [fn, index] = getVariableField@PhysicalModel(model, name);
                varfound = true;
            end

        end

        %%
        function fds = getAllVarsNames(model)
            fds = {'components', 'masterComponents', 'logcomponents', ...
                   'logmasterComponents', model.CompNames{:}, ...
                   model.MasterCompNames{:}, model.logCompNames{:}, ...
                   model.logMasterCompNames{:}};
        end

        %%
        function [state, report, model] = initState(model, userInput, varargin)
        %initState solves the chemical system using the specified input and
        % information inside of the already defined state.
        %
        % SYNOPSIS:
        %  [state, report, model] = initState(model, userInput, varargin)
        %
        % DESCRIPTION:
        %   A function used to solve the chemical system generated by
        %   ChemModel according to constraints given in userInput. 
        %
        % REQUIRED PARAMETERS:
        %   userInput        - A matrix containing the values of elements
        %       or species specified as inputs during the generation of the
        %       chemical system using ChemModel. userInput must have as
        %       mnay columns as there are specified knowns, which is equal
        %       to the number of entries in elements. 
        %           
        %
        % OPTIONAL PARAMETERS:
        %   state            - A structure containg information about the
        %       physical model. Currently the only relevant information for 
        %       the chemical system is state.temperature. This is an
        %       optional parameter, when not provided temperature is
        %       assumed to be 298 K. Should be given as a 'key'/value pair.
        %
        %           state.temperature = 170*Kelvin;
        %           
        %           state = model.initState(userInput, 'state', state)
        %
        %   chargeBalance   - This input toggles the enforcment of strict
        %       charge balance. A master component must be specified. If
        %       charge balance is not satisfied with the given constraint
        %       the solver will adjust the total concentration of the
        %       specified element to enforce charge balance. chargeBalance
        %       must be given as a 'key'/value pair. The specified element
        %       must be one of the elements provided as an input. By
        %       defualt charge balance is not strictly enforced. This make
        %       the chemical solver more robust. 
        %
        %           state = model.initState(userInput, 'chargeBalance', 'Cl')
        %
        %
        % OUTPUTS:
        %   state           - A structure containing the solution to the
        %       chemical system. This includes state.masterCompoents and
        %       state.components. masterComponents contains the total
        %       element concentrations and surface functional group
        %       concentrations. components contains the species
        %       concentrations as well as the surface acitivity
        %       multipliers. Values can be retreived from state using the
        %       model.getProp() function.
        %
        %   report          - The nonlinear solver report for the chemical
        %       system. 
        %
        %   model           - The updated chemical model object.
        %
        % SEE ALSO:
        %   ChemicalModel
            
            
            p = inputParser;
            
            valInd = cellfun(@(x) isempty(x), regexpi(model.MasterCompNames, '>'));
            
            valFun = @(x) any(validatestring(x, model.MasterCompNames(valInd)));
            p.addOptional('chargeBalance', 'nochargebalance', valFun);
            p.addOptional('state',struct,@isstruct)
            
            p.parse(varargin{:})
            
            givenTest = any(strcmpi(p.Results.chargeBalance, horzcat(model.chemicalInputModel.inputNames,'nochargebalance')));
            assert(givenTest, ['Only elements whos values are given (marked with "*") can be used for charge balance.']);
            
            chargeBalance = ~strcmpi(p.Results.chargeBalance,'nochargebalance');

            model.chemicalInputModel = model.chemicalInputModel.validateModel();
            model.chemicalInputModel.plotIter = model.plotIter;
            
            if isempty(model.compositionModel)
                model.compositionModel = compositionModel();
            end
            if isempty(model.compositionReactionModel)
                model.compositionReactionModel = compositionReactionModel();
            end
            props = properties(model.chemicalInputModel);
            for i = 1 : numel(props);
            	model.compositionModel.(props{i}) = model.chemicalInputModel.(props{i});
                model.compositionReactionModel.(props{i}) = model.chemicalInputModel.(props{i});
            end
                
            
            model.compositionModel = model.compositionModel.validateModel();
            model.compositionReactionModel = model.compositionReactionModel.validateModel();

            cheminput = model.chemicalInputModel;
            inputNames = cheminput.inputNames;
            unknownNames = cheminput.unknownNames;

            inSize = size(userInput);
            if ~isempty(model.surfInfo);
                k = numel(model.surfInfo.master);
            else
                k = 0;
            end
            
            assert(inSize(2) == model.nMC-k, ['For the specified chemical system the input constraint must have ' ...
                                num2str(model.nMC-k) ' columns.']);

            
            state.masterComponents = eps*ones(size(userInput,1), model.nMC);
            state.components = eps*ones(size(userInput,1), model.nC);
            state.combinationComponents = eps*ones(size(userInput,1), model.nLC);
            
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
            
            % make sure they are the correct size
            givenSize = cellfun(@(x) size(x,1), in);
            mSize = max(givenSize);
            
            assert(any(givenSize == mSize | givenSize == 1), 'Input size must be consistent between input and surface parameters');
           
            for i = 1 : model.nMC
            	state = model.setProp(state, inputNames{i}, in{i}.*ones(mSize,1));
            end
            
            state = model.syncLog(state);
            
            % create initial guess
            fprintf('Computing initial guess...\n')
            [state, ~, report_c] = model.compositionModel.solveChemicalState(state);
            [state, ~, report_cr] = model.compositionReactionModel.solveChemicalState(state);
            
            % solve chemical system
            fprintf('Solving chemical system...\n')
            [state, ~, report] = model.chemicalInputModel.solveChemicalState(state);
            
            state0 = state;
            
            if chargeBalance
                fprintf('Enforcing strict charge balance...\n');
                if isempty(model.chargeBalanceModel)
                    model.chargeBalanceModel = chargeBalanceModel();
                end
                for i = 1 : numel(props);
                    model.chargeBalanceModel.(props{i}) = model.chemicalInputModel.(props{i});
                end
                model.chargeBalanceModel.CVC = p.Results.chargeBalance;
                [state, ~, report] = model.chargeBalanceModel.solveChemicalState(state);
                
                if size(state.components,2) ~= numel(model.CompNames)
                    warning('Convergence not acheived for charge balance with given inputs. Try a different element for charge compensation, or increase salt content.');
                    state = state0;
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
            % find mass input components and remove the asterik from the name
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

            % create the index matricies for construction of secondary components
            tmpvec = zeros(nM, 1);
            for i = 1 : nM
                model.MCind{i} = tmpvec;
                model.MCind{i}(i) = 1;
            end

            % make a list of the component names
            model.MasterCompNames = cell(1,nM);
            [model.MasterCompNames{:}] = deal(names{:});

            model.logMasterCompNames = cellfun(@(name) ['log', name], ...
                                               model.MasterCompNames, ...
                                               'uniformoutput', false);

            % write number of master components
            model.nMC = nM;

        end

        %%
        function [model] = initSecondaryComponents(model, secondarycomponents )
            %initSecondarySpecies computes secondary components as linear combinations
            %of master components
         

            assert(~isempty(model.MasterCompNames),'The function initSecondaryComponents can only be called after initMasterComponents.');

            % make sure input is a name value pair of the correct size
            assert(iscell(secondarycomponents),'input must be a cell array of strings');

            names = regexprep(secondarycomponents,'[\s]','','ignorecase');

            test = cellfun(@(x) ~isempty(x), regexpi(secondarycomponents, '(psi)(sig)'));
            assert(sum(test) == 0, 'The use of "sig" and "psi" are reserved for use within ChemicalModel, please remove them from the species cell');
            
            % find input components and remove the asterik from the name
            InInd = cellfun(@(x) ~isempty(x), regexp(names, '*'));
            names = regexprep(names, '*','');



            if isempty(model.chemicalInputModel)
                model.ChemicalInputModel = ChemicalInputModel();
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
                assert(~(sum(strcmpi(names{i}, model.MasterCompNames)) > 0), ['The chemical species ' '"' names{i} '" is a repeat of an element or surface functional group.']);
            end

            % length of master component string
            strlen = cellfun(@length,model.MasterCompNames);

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
            model.ChargeVector = tmpvec;
            surfInd = cellfun(@(x) ~isempty(x), regexpi(namesO,'>'));
            model.ChargeVector(surfInd) = 0;

            % replace master component names with index vectors
            for i = 1 : nS;
                names{i} = regexprep(names{i}, '[-+](\d*([./]\d*)?)*','','ignorecase');
                for j = 1 : nM
                    ind = I(j);
                    [match,nomatch] = regexpi(names{i}, '(model\.MCind{\d})','match','split');
                    nomatch = strrep(lower(nomatch), lower(model.MasterCompNames{ind}), ['model.MCind{' num2str(ind) '}']);
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
            end

            tmpvec = zeros(nS, 1);

            for i = 1 : nS
                model.Cind{i} = tmpvec;
                model.Cind{i}(i) = 1;
            end

            % store cell array of secondary component names
            model.CompNames = cell(1,nS);
            [model.CompNames{:}] = deal(namesO{:});

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
            model.CompositionMatrix = horzcat(indsum{:});

            model.logCompNames = cellfun(@(name) ['log', name], model.CompNames, ...
                                         'uniformoutput', false);

            model.nC = numel(model.CompNames);
            
            for i = 1 : numel(model.MasterCompNames)
                MCName = model.MasterCompNames{i};
                ind = cellfun(@(x) ~isempty(x), regexpi(model.CompNames,MCName));
                assert(sum(ind) > 0, ['The element ' MCName ' appears to make no contribution to a chemical species. Verify species and element names and/or remove the element if it is not needed.']); 
            end
            


        end

        %%
        function [model] = initChemicalReactions(model, varargin)
            %initChemicalReactions computes rection matrix as linear combinations of
            %secondary components, must be called after initMasterComponents and
            %initSecondaryComponents
           
            assert(~isempty(model.CompNames),'The function initChemicalReactions can only be called after initMasterComponents and initSecondaryComponents');

            % make sure input is a name value pair
            assert(rem(numel(varargin),2) == 0, 'Reaction and reaction constant must be given as a name/value pair.');

            rxn = regexprep(varargin(1:2:end),'[\s]','');

            rxnK = varargin(2:2:end);

            % unwrap master component vectors
            nRx = numel(rxn);
            nS = model.nC;
            nM = model.nMC;

            % check that an appropriate number of inputs have been
            % designated
            strlen = cellfun(@length,model.CompNames);
            [~, I] = sort(strlen, 2, 'descend');

            % make a copy of reactions
            rxnsO= rxn;

            % rewrite reactions in terms of the master component index
            % vectors
            for i = 1 : nRx;
                for j = 1 : nS
                    ind = I(j);
                    [match,nomatch] = regexpi(rxn{i}, '(model\.Cind{\d})','match','split');
                    nomatch = strrep(lower(nomatch), lower(model.CompNames{ind}), ['model.Cind{' num2str(ind) '}']);
                    rxn{i} = strjoin(nomatch,match);
                end
                rxn{i} = strrep(rxn{i}, 'cind', 'Cind');
                remain = regexpi(rxn{i}, '[^(model\.Cind{\d})+-/*\.\d(<->)]','once');
                assert(isempty(remain),['The chemical reaction "' rxnsO{i} '" appears to contain a string combination that does not correspond to secondary components. Ensure the species names contained in the reaction are consistent with those listed in secondarycomponents.']);
            end

            % check for <-> and handle it appropriately
            for i = 1 : nRx
                assert(~isempty(regexpi(rxn{i},'(<->)','once')),...
                    ['Chemical reaction ' num2str(2) ' is improperly formatted. Reactions must be written as "reactants <-> products". If no reactant or product exists place a "0" on the appropirate side.']);
                rxn{i} = regexprep(rxn{i}, '<->', ')+','ignorecase');
                try
                    tmp{i} = eval(['-(' rxn{i} ,';'])';
                catch
                    error(['The reaction "' rxnsO{i} '" appears to contain a string combination that does not correspond to species names. Make sure to add all relevant chemical species to the species list. If a reaction involves more than one species use "*" to designate this (i.e. 2*Na+ rather than 2Na+). ']);
                end
            end

            % create reaction matrix
            model.ReactionMatrix = vertcat(tmp{:});

            % assemble reaction constants
            model.ReactionConstants = [rxnK{:}]';
            model.LogReactionConstants = log(model.ReactionConstants);

            % check that there are an equal number of each master component on both sides of the chemical reacitons.
            for i = 1 : nRx
                mCind = model.ReactionMatrix(i,:) ~= 0;
                mCcomp = repmat(model.ReactionMatrix(i,mCind),nM,1).*model.CompositionMatrix(:,mCind);
                balance = sum(sum(mCcomp,2));
                assert(balance == 0, ['The chemical reaction "' rxnsO{i} '" does not balance elements and/or surface functional groups.'])
            end
            
            % check for charge balance within each reaction
            for i = 1 : nRx
                mCind = model.ReactionMatrix(i,:) ~= 0;
                mCcomp = model.ReactionMatrix(i,mCind).*model.allCharge(:,mCind);
                balance = sum(sum(mCcomp,2));
                assert(balance == 0, ['The chemical reaction "' rxnsO{i} '" does not balance charge.'])
            end
            
            model.nR = nRx;
            model.rxns = rxnsO;

            surfInfo = model.surfInfo;


            RM = model.ReactionMatrix;

            nP = sum(cellfun(@(x) ~isempty(x), regexpi(model.CompNames,'psi')));
            RM = horzcat(RM, zeros(model.nR, nP));

            % Add surface activities to reaction matrix and remove
            % contribution to aqueous charge conservation
            if model.surfFlag && ~isempty(model.surfInfo)
                
                gNames = model.surfaces.groupNames;
                mNames = model.surfaces.groupMasterNames;
                
                for i = 1 : numel(gNames)
                    for j = 1 : numel(mNames(i,:));
                        iterNames = mNames(cellfun(@(x) ~isempty(x), mNames(i,:)));
                        iterInd = strcmpi(iterNames{j}, model.surfInfo.master);
                        sNames = surfInfo.species{iterInd};
                        switch model.surfaces.scm{i}
                            case {'tlm','ccm'}
                                % find the number of layers
                                nL = numel(surfInfo.charge{iterInd}{1});
                                layerInd = cellfun(@(x) ~isempty(x), regexpi(model.CompNames,[gNames{i} '_']));
                                for k = 1 : numel(sNames)
                                    % find the columns of the reaction matrix that contain
                                    % the species
                                    sp = sNames{k};
                                    spInd = strcmpi(sp, model.CompNames);
                                    rmVec = RM(:,spInd);

                                    % replace the correct column and row with the
                                    % contribution
                                    RM(logical(rmVec),layerInd) = RM(logical(rmVec),layerInd) + repmat(rmVec(rmVec ~= 0),1,nL).*repmat(surfInfo.charge{iterInd}{k},sum(logical(rmVec)),1);
                                end
                        end
                    end
                end
            end

            
            nPsi = sum(cellfun(@(x) ~isempty(x), regexpi(model.CompNames, 'psi')));
            assert(model.nR + model.nMC == model.nC - nPsi, ['The given chemical system is rank deficient. The number of species must be equal to the sum of the number elements, surface functional groups, and reactions.'])
            
            nI = numel(model.chemicalInputModel.inputNames);
            
            assert(nI == model.nMC , ['For the defined chemical system ' num2str(model.nMC) ' components, elements or species, must be designated as inputs. Use an "*" to designated a component as an input.']);
            model.ReactionMatrix = RM;


        end

        function model = setupMaxMatrices(model)

            C = model.CompositionMatrix;
            nC = model.nC;
            maxMatrices = cell(nC, 1);
            for i = 1 : nC
                maxMatrices{i} = diag(1./C(:, i));
            end
            model.maxMatrices = maxMatrices;

        end


        %%
        function model = initLinearCombinations(model, linearCombos)
            %initElectrostaticModel parses the surface chemistry inputs 
            
            newVar = linearCombos(1:2:end);
            combos = linearCombos(2:2:end);
            
            nVar = numel(newVar);
            
            % find input components and remove the asterik from the name
            InInd = cellfun(@(x) ~isempty(x), regexp(newVar, '*'));
            names = regexprep(newVar, '*','');



            if isempty(model.chemicalInputModel)
                model.ChemicalInputModel = ChemicalInputModel();
            end
            model.chemicalInputModel.inputNames = ...
                horzcat(model.chemicalInputModel.inputNames, names(InInd));
            
            
            % unwrap master component vectors
            nS = model.nC;
            nM = model.nMC;

            % sort comp names by length
            strlen = cellfun(@length,model.CompNames);
            [~, I] = sort(strlen, 2, 'descend');
            
            combos0 = combos;
            
            % rewrite reactions in terms of the master component index
            % vectors
            for i = 1 : nVar;
                for j = 1 : nS
                    ind = I(j);
                    [match,nomatch] = regexpi(combos{i}, '(model\.Cind{\d})','match','split');
                    nomatch = strrep(lower(nomatch), lower(model.CompNames{ind}), ['model.Cind{' num2str(ind) '}']);
                    combos{i} = strjoin(nomatch,match);
                end
                combos{i} = strrep(combos{i}, 'cind', 'Cind');
            end 
            
            % check for <-> and handle it appropriately
            for i = 1 : nVar
                try
                    tmp{i} = eval([ combos{i} ,';'])';
                catch
                    error(['The linear combination "' combos0{i} '" appears ' ...
                                        'to contain a string combination that ' ...
                                        'does not correspond to species names. ' ...
                                        'Make sure to add all relevant chemical ' ...
                                        'species to the species list. Use  ' ...
                                        '"*" after the coefficients in the ' ...
                                        'linear combination  (i.e. 2*Na+ rather than 2Na+). ']);
                end
            end
            
            % create reaction matrix
            model.CombinationMatrix = vertcat(tmp{:});
            model.CombinationNames = names;
            model.nLC = numel(names);
            
           

        end
%%
        function model = initElectrostaticModel(model, givenSurfaceInformation, varargin)
            %initElectrostaticModel parses the surface chemistry inputs 
            
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
                
                % make sure given species are correct
                Mind = strcmpi(MName,model.MasterCompNames);
                ind = logical(model.CompositionMatrix(Mind,:));
                speciesNames{i} = model.CompNames(ind);
                
                
                % grab the charge from the species name
                nS = numel(speciesNames{i});
                tmpvec = zeros(1, nS);
                
                for j = 1 : nS
                    sInd = strcmpi(speciesNames{i}{j}, model.CompNames);
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
                n = 0;
                
                mastNames = model.surfaces.groupMasterNames(i,:);
                
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
                assert(isequal(mTest{:},mTest{:}), 'Surface functional groups that are combined into a single electrostatic surface must use the same electrostatic model.');
                assert(isequal(cTest{:},cTest{:}), 'Surface functional groups that are combined into a single electrostatic surface must use the same values for capacitance.');
                
                model.surfaces.scm{i} = mTest{1};
                model.surfaces.c{i} = cTest{1};                

                switch mTest{1}
                    case 'ccm'
                        model.CompNames{end+1} = [ m '_ePsi'];
                        model.surfaceChargeNames{end+1} = [m '_sig'];
                        model.surfacePotentialNames{end+1} = [m '_Psi'];
                        n = 1;
                    case 'tlm'
                        model.CompNames =  horzcat( model.CompNames, [m '_ePsi_0'], [m '_ePsi_1'], [m '_ePsi_2']);
                        model.surfaceChargeNames= horzcat(model.surfaceChargeNames, [m '_sig_0'],[m '_sig_1'],[m '_sig_2']);
                        model.surfacePotentialNames = horzcat(model.surfacePotentialNames, [m '_Psi_0'],[m '_Psi_1'],[m '_Psi_2']);
                        n = 3;
                end
                model.CompositionMatrix = [model.CompositionMatrix, zeros(model.nMC, n)];
            end

            model.logCompNames = cellfun(@(name) ['log', name], model.CompNames, ...
                                         'uniformoutput', false);

            model.nC = numel(model.CompNames);


        end
        
        
        %%
        function printChemicalSystem(model, varargin)
            %printchemicalSystem outputs the chemical system tableua
            
            chemicalSystemPrintFunction(model, varargin{:});
          

        end


        %%
        function [state, report] = updateState(model, state, problem, dx, drivingForces) %#ok
        % Update state based on Newton increments
            [state, report] = updateState@PhysicalModel(model, state, problem, ...
                                                        dx, drivingForces);

            state = model.syncFromLog(state);
            
            surfParam = sum(cellfun(@(x) ~isempty(x) , regexpi(model.CompNames, 'psi'))); 
                        
            
            nonLogVariables = regexprep(problem.primaryVariables, 'log', '');

            if surfParam > 0
            	[names, mins, maxs] = computeMaxPotential(model, state); 
            end
            
            len = cellfun(@(x) length(x), nonLogVariables);
            [~,sortInd] = sort(len(:),1, 'ascend');
            pVar = nonLogVariables(sortInd);
            
            for i = 1 : numel(pVar)
                
                p = pVar{i};
                compInd = strcmpi(p, model.CompNames);
                
                if any(strcmpi(p, model.MasterCompNames))
                    state = model.capProperty(state, p, 1e-50, 2.5*mol/litre);
                elseif ~isempty(regexpi(p, 'psi'))
                	ind = find(strcmpi(p, names));
                    state = model.capProperty(state, p, mins{ind}, maxs{ind});
                elseif ismember(p, model.CombinationNames)
                    state = model.capProperty(state, p, -2.5*mol/litre, 2.5*mol/litre);
                else
                    maxvals = model.maxMatrices{compInd}*((state.masterComponents)');
                    maxvals = (min(maxvals))';             
                    state = model.capProperty(state, p, 1e-50, maxvals); 
                end
                
            end
            
            state = model.syncLog(state);
            

            if model.plotIter 
                h = findobj('tag', 'updatechemfig');
                if isempty(h)
                    figure
                    set(gcf, 'tag', 'updatechemfig');
                    h = findobj('tag', 'updatechemfig');
                end
                set(0, 'currentfigure', h)
                clf
                if size(state.components, 1) ~= 1
                    plot(log10(state.components*litre/mol));
                    title('components (chemistry step)');
                    legend(model.CompNames);
                end

                h = findobj('tag', 'chemmastercomp');
                if isempty(h)
                    figure
                    set(gcf, 'tag', 'chemmastercomp');
                    h = findobj('tag', 'chemmastercomp');
                end
                set(0, 'currentfigure', h)
                clf
                if size(state.components, 1) ~= 1
                    plot(log10(state.masterComponents*litre/mol));
                    title('master components (chemistry step)');
                    legend(model.MasterCompNames);
                end
                drawnow;
            end

        end

        %%
        function state = validateState(model, state)
            state = validateState@PhysicalModel(model, state);
            if (~isfield(state, 'logcomponents') || ~isfield(state, 'logmastercomponents') )
                state = model.syncLog(state);
            end
            state = model.syncLog(state);
        end

        %%
        function state = syncLog(model, state)
            state.logmasterComponents = log(state.masterComponents);
            state.logcomponents       = log(state.components);
        end

        %%
        function state = syncFromLog(model, state)
            state.masterComponents = exp(state.logmasterComponents);
            state.components       = exp(state.logcomponents);
        end


        %%
        function [logComps, logMasterComps] = prepStateForEquations(model, state)

            logComps = cell(model.nC, 1);
            [logComps{:}] = model.getProps(state, model.logCompNames{:});

            logMasterComps = cell(model.nMC, 1);
            [logMasterComps{:}] = model.getProps(state, model.logMasterCompNames{:});

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
        %computeAcitivities computes the acitivity of each aqueous species
        %using the extended Davies equaiton. 
        %
        % SYNOPSIS:
        %  [state] = computeActivities(model, state)
        %
        %
        % REQUIRED PARAMETERS:
        %   state        - the state variable produced by model.initState.
        %       Must at least include the field state.components. 
        %          
        %
        % OUTPUTS:
        %   state           - A structure containing the acitivity of each
        %   aqueous species in units of mol/meter^3. Activities can be
        %   retrieved using the getProp command by calling for the species
        %   name prepended by 'a.'
        %
        % EXAMPLE:
        %
        %   state = chem.computeSurfaceCharges(state);
        %   aH2O = chem.getProp(state, 'aH2O');
        %
        %
        % SEE ALSO:
        %   computeChargeBalance
        
            assert( nargout == 2, ['Output argument to computeActivities must include the state variable and the chemical model object']);
            [state, model] = activity(model, state);
            
        end
        
        %%
        function [state, model] = computeChargeBalance(model, state)
        %computeChargeBalance computes the residual of the aqueous charge
        %balance equation.
        %
        % SYNOPSIS:
        %  [state] = computeChargeBalance(model, state)
        %
        %
        % REQUIRED PARAMETERS:
        %   state        - the state variable produced by model.initState.
        %       Must at least include the field state.components. 
        %          
        %
        % OUTPUTS:
        %   state           - A structure containing the value of the
        %   aqueous charge balance residual in units of percent of total
        %   charge. The charge balance can be retrieved using the getProps
        %   command by calling for the variables 'chargeBalance'.
        %
        %
        %   chem            - updated chemical object now includes the names 
        %       for calling the charge balance using getProps 
        %
        % EXAMPLE:
        %
        %   state = chem.computeSurfaceCharges(state);
        %   charge = chem.getProp(state, 'chargeBalance');
        %
        %
        % SEE ALSO:
        %   computeActivities

            assert( nargout == 2, ['Output argument to computeChargeBalance must include the state variable and the chemical model object']);
            [state, model] = chargeBalance(model, state);
            
        end
        
        %%
        function [state, model] = computeSurfacePotentials(model, state)
        %computeSurfacePotentials computes the potential of each layer of each
        %surface and adds the values to the field state.surfacePotentials.
        %
        % SYNOPSIS:
        %  [state] = computeSurfacePotentials(model, state)
        %
        %
        % REQUIRED PARAMETERS:
        %   state        - the state variable produced by model.initState.
        %       Must at least include the field state.components. 
        %          
        %
        % OUTPUTS:
        %   state           - A structure containing the potential of each layer of
        %       each surface in Volts. Values can be retrieved using
        %       the getProps command by calling for the surface functional
        %       group name, followed by '_Psi_' followed by the layer
        %       number (0, 1, 2) or for the constant capacitance model just
        %       use '_Psi'.
        %
        %
        %   chem            - updated chemical object now includes the names 
        %       for calling the surface potentials using getProps 
        %
        % EXAMPLE:
        %
        %   state = chem.computeSurfaceCharges(state);
        %   potential0 = chem.getProp(state, '>SiO_Psi_0');
        %   potential1 = chem.getProp(state, '>SiO_Psi_1');
        %
        %
        % SEE ALSO:
        %   computeSurfaceCharges
        
            assert( nargout == 2, ['Output argument to computeSurfacePotentials must include the state variable and the chemical model object']);
            [state, model] = surfacePotential(model, state);
            
        end
        
        %%
        function [state, model] = computeSurfaceCharges(model, state)
        %computeSurfaceCharges computes the charge of each layer of each
        %surface and adds the values to the field state.surfaceCharges.
        %
        % SYNOPSIS:
        %  [state] = computeSurfaceCharges(model, state)
        %
        %
        % REQUIRED PARAMETERS:
        %   state        - the state variable produced by model.initState.
        %       Must at least include the field state.components. 
        %          
        %
        % OUTPUTS:
        %   state           - A structure containing the charge of each laye of
        %       each surface in C/meter^2. Values can be retrieved using
        %       the getProps command by calling for the surface functional
        %       group name, followed by '_sig_' followed by the layer
        %       number (0, 1, 2) or for the constant capacitance model just
        %       use '_sig'.
        %
        %
        %   chem            - updated chemical object now includes the names 
        %       for calling the surface charges using getProps 
        %
        % EXAMPLE:
        %
        %   state = chem.computeSurfaceCharges(state);
        %   charge0 = chem.getProp(state, '>SiO_sig_0');
        %   charge1 = chem.getProp(state, '>SiO_sig_1');
        %
        %
        % SEE ALSO:
        %   computeSurfacePotentials
        
            assert( nargout == 2, ['Output argument to computeSurfaceCharges must include the state variable and the chemical model object']);                
            [state, model] = surfaceCharge(model, state);
        end
        
        %%
        function [state, model] = computeAqueousConcentrations(model, state)
        %computeAqueousConcentrations computes the total concentration of each 
        % element in the aqueous phase.
        %
        % SYNOPSIS:
        %  [state] = computeAqueousConcentrations(model, state)
        %
        %
        % REQUIRED PARAMETERS:
        %   state        - the state variable produced by model.initState.
        %       Must at least include the field state.components. 
        %          
        %
        % OUTPUTS:
        %   state           - A structure containing the total aqueous
        %   concentration of each element in the field
        %   aqueousConcentrations. The values can be retrieved using the
        %   getProps function by asking for the element name followed
        %   by '_aq'.
        %
        %   chem            - updated chemical object now includes the names 
        %       for calling the aqueous concetrations using getProps 
        %
        % EXAMPLE:
        %
        %   state = chem.computeAqueousConcentrations(state);
        %   Na = chem.getProp(state, 'Na_aq');
        %   H = chem.getProp(state, 'H_aq');
        %
        %
        % SEE ALSO:
        %   computeSurfaceConcentrations
        
            assert( nargout == 2, ['Output argument to computeAqueousConcentrations must include the state variable and the chemical model object']);
            [state, model] = aqueousConcentrations(model, state);
            
        end
        
        %%
        function [state, model] = computeSurfaceConcentrations(model, state)
        %computeSurfaceConcentrations computes the total concentration of each 
        % element on the surface.
        %
        % SYNOPSIS:
        %  [state] = computeAqueousConcentrations(model, state)
        %
        %
        % REQUIRED PARAMETERS:
        %   state        - the state variable produced by model.initState.
        %       Must at least include the field state.components. 
        %          
        %
        % OUTPUTS:
        %   state           - A structure containing the total surface
        %       concentration of each element in the field
        %       surfaceConcentrations. The values can be retrieved using the
        %       getProps function by asking for the element name followed
        %       by '_surf'.
        %
        %   chem            - updated chemical object now includes the names 
        %       for calling the surface concetrations using getProps 
        %
        % EXAMPLE:
        %
        %   [state chem] = chem.computeAqueousConcentrations(state);
        %   Na = chem.getProp(state, 'Na_surf');
        %   H = chem.getProp(state, 'H_aq');
        %
        %
        % SEE ALSO:
        %   computeSurfaceConcentrations
        
            assert( nargout == 2, ['Output argument to computeSurfaceConcentrations must include the state variable and the chemical model object']);
            [state, model] = surfaceConcentrations(model, state);
            
        end
    end

end
