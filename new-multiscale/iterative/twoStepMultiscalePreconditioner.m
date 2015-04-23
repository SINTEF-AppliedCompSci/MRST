function x = twoStepMultiscalePreconditioner(A, b, solveCoarse, smooth, it)
    if nargin < 5
        it = 1;
    end
    if it == 1
        x = smooth(b);
        x = x + solveCoarse(b - A*x);
    else
        x = zeros(size(b));
        for i = 1:it
            x = x + smooth(b - A*x);
            x = x + solveCoarse(b - A*x);
        end
    end
end
