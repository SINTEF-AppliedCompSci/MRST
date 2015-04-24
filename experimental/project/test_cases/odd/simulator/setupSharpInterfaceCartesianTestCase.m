function testCase = setupSharpInterfaceCartesianTestCase(CO2, varargin)

    opt = testCaseDefaults();                       % load non-grid related defaults 
    opt = concatStruct(opt, cartesianGridDefaults); % include grid-related defaults
    opt = merge_options(opt, varargin{:});          % override with user options

    % Ensure any scalar parameters are expanded into vectors of the
    % appropriate size, and verify results.
    opt.topo       = expand_var(opt.topo,     opt.cartDims(1:2)+1);
    opt.porosity   = expand_var(opt.porosity, opt.cartDims);
    opt.perm       = expand_var(opt.perm,     prod(opt.cartDims));

    %% Preparing grid structure
    G = cartGrid(opt.cartDims, opt.physDims);
    
    % adjusting topography 
    G.nodes.coords(:,3) = ...
        G.nodes.coords(:,3) + repmat(opt.topo(:), G.cartDims(3)+1, 1);

    % computing geometry
    G = computeGeometry(G);
    
    %% Preparing rock structure
    rock.perm = opt.perm;
    rock.poro = opt.porosity;
    
    %% Remove fields of 'opt' that 'setupSharpInterfaceTestCase' does
    %  not expect.
    opt = rmfield(opt, 'cartDims');
    opt = rmfield(opt, 'physDims');
    opt = rmfield(opt, 'topo');
    opt = rmfield(opt, 'perm');
    opt = rmfield(opt, 'porosity');

    %% Calling the function in charge of setting up all non-grid parameters
    cell_opt = fullStructToCell(opt);
    testCase = setupSharpInterfaceTestCase(CO2, G, rock, cell_opt{:});
    
end 
            
            