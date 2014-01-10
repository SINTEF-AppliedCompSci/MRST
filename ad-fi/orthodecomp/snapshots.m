function [X X_bar] = snapshots(values)
    X_bar = mean(values,2);
    X = values - repmat(X_bar, 1, size(values,2));
end
