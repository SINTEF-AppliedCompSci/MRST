function [isOscillating, isStagnating] = detectNewtonOscillations(history, primary, current, tol)
    if current < 3
        isOscillating = false;
        isStagnating = false;
        return
    end

    tmp = history(current-2:current, primary);

    oscillate =  relChange(tmp(1,:), tmp(3,:)) < tol & ...
                 relChange(tmp(2,:), tmp(3,:)) > tol;

    stagnate = relChange(tmp(2,:), tmp(1,:));
    stagnate(isnan(stagnate)) = 0;

    isStagnating = all(stagnate < 1e-3);
    isOscillating = sum(oscillate) > 1;
end

function v = relChange(a,b)
    v = abs((a-b)./b);
end
