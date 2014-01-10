function p = msfvPressureSolve(arg1, rhs, DG, useCorrection)
    persistent B Cx updateTargets clusters A_nn;
    if nargin == 0
        % Reset all basis functions
        B = [];
        A_nn = [];
        Cx = [];
        return
    end

    if nargin == 1
        % Flag cells for update before next solve
        updateTargets = arg1;
        return
    end

    A = arg1;
    nc = size(A, 1);

    if isempty(B)
        if useCorrection && isempty(Cx)
            [B, Cx] = createMSFVBasis(A, DG, useCorrection);
        else
            B = createMSFVBasis(A, DG, useCorrection);
        end
        updateTargets = false(nc, 1);
    end

%     if any(updateTargets)
%         % Do update
%         if isempty(clusters)
%
%             clusters = zeros(nc, 1);
%             ni = numel(DG.ii);
%
%             Ap = DG.P*A*DG.P';
%             Ap = Ap(1:ni, 1:ni);
%             clusters(DG.ii) = components(Ap);
%         end
%         B = createMSFVBasis(A, DG, updateTargets, clusters);
%         updateTargets = false(size(A, 1), 1);
%     end
    X = DG.X;

    if isempty(A_nn)
        A_nn = X*A*B;
    end
    if useCorrection
        C_rhs = Cx(rhs);
        q = X*(rhs - A*C_rhs);
    else
        C_rhs = zeros(size(rhs));
    end

    p = B*mldivide(A_nn, q) + C_rhs;

end
