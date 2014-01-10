function printResidual(residuals, gmresits, eqnnames, iteration, CNV, MB)
    if iteration == 1
        fprintf('%-9s', eqnnames{:})
        fprintf('\n');
    end
    CNV = CNV(CNV ~= 0);
    MB  = MB(CNV ~= 0);
    fprintf('%8.2e ', residuals);
    fprintf('** CNV: ');
    fprintf('%2.2e ', CNV);
    fprintf('MB: ');
    fprintf('%2.2e ', MB);
    if ~isempty(gmresits)
        fprintf('** %d GMRES iterations', gmresits(2));
    end
    fprintf('\n');
end
