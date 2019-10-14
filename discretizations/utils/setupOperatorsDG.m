function op = setupOperatorsDG(disc, G, rock, varargin)

    op = setupOperatorsTPFA(G, rock, varargin{:});
    
    [~, ~, ~, f] = disc.getCubature(find(disc.internalConn), 'face');
    nf = numel(f);
    
    N = [1:nf; nf+1:2*nf]';
    
    M = sparse((1:nf)'.*[1,1], N, 0.5, nf, 2*nf);
    op.M = M;
    op.faceAvg = @(x) M*x;
    
    op.faceUpstr = @(flag, x) faceUpstr(flag, x, N, [nf, 2*nf]);
    
    C = sparse((1:nf)'.*[1,1], N, ones(nf,1)*[1 -1], nf, 2*nf);
    op.C = C;
    op.Grad = @(x) -C*x;
    
    op.Div = @(x) -C*x;
    
    op.T = op.T_all(f);
    
    op.velocityInterp = velocityInterpolation(G, 'mimetic');
    
end

