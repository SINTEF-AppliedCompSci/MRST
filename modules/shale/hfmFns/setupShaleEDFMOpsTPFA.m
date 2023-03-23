function s = setupShaleEDFMOpsTPFA(G, rock, tol, varargin)
    % Modifies from setupEDFMOperatorsTPFA.m.
    
    opt = struct('deck', [], 'neighbors', [], 'trans', [], 'porv', []);
    opt = merge_options(opt, varargin{:});

    s = setupEDFMOperatorsTPFA(G, G.rock, tol);

    T = getFaceTransmissibility(G, rock, opt.deck);
    s.T_all = T;
    T= [T; G.nnc.T];
    T = T(s.internalConn);
    s.T = T;

    % HR-EDIT- Include adsorption
    if isfield(G.rock, 'shaleMechanisms') && isfield(G.rock.shaleMechanisms,'sorption')
       s.isSorbed = rock.isMatrix;
       s.gv = s.isSorbed.*(1-rock.poro) .* G.cells.volumes;
    end
    % Identification of matrix with value 1
    %OMO edit: test that isMatrix exists
    if isfield(rock, 'isMatrix')
        s.isMatrix=rock.isMatrix;
    end

    %OMO edit: migrate hfm from 2018b to 2019b
    n = size(s.N,1);
    C = (s.C)';  % transpose it once for speed
    if n == 0	
        % We have zero faces. Account for Matlab's preference for	
        % reducing expressions of type a + [] to [].	
        s.AccDiv = @(acc, flux) acc;	
        s.Div = @(x) zeros(G.cells.num, 1);	%nc = G.cells.num = max(max(N))
    else	
        s.AccDiv = @(acc, flux) acc + C*flux;	
        s.Div = @(x) C*x;	
    end

end

