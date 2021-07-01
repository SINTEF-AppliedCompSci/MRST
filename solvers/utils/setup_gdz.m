function gdz = setup_gdz(model, grad_op)

    g = model.getGravityVector();
    griddim = model.G.griddim;
    grad = zeros(sum(model.operators.internalConn), griddim);

    % FIX move loop
    for k = 1:griddim
        grad(:, k) = grad_op(model.G.cells.centroids(:, k));
    end

    gdz = grad * g';
end