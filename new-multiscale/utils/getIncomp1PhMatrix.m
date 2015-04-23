function A = getIncomp1PhMatrix(G, T, state, fluid)
    
    if nargin < 3
        state = initResSol(G, 0);
    end
    
    if nargin < 4
        fluid = initSingleFluid('rho', 1, 'mu', 1);
    end
    
    use_trans = numel(T) == G.faces.num;
    state = incompTPFA(state, G, T, fluid, 'MatrixOutput', true, ...
                    'LinSolve', @(A, x) 0*x, 'use_trans', use_trans);
    A = state.A;
    % Undo magic scaling from the inside of incompTPFA
    A(1, 1) = A(1, 1)/2;
end
