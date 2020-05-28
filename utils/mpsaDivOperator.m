function divu = mpsaDivOperator(sol, extforce, R1, R2, div)
    % Compute nodeface displacement 
    unf = R1*sol + R2*extforce;
    % Compute divergence
    divu = div*unf;
end
