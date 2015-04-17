function u=fluidDerivative(fluid,p,T,property)
% calculate joule tomson coefficent
[pa,Ta] = initVariablesADI(double(p),double(T));
switch property
    case 'jouleThomson'
        val=fluid.enthalpy(pa,Ta);
        u=-full(diag(val.jac{1}))./full(diag(val.jac{2}));
    case 'freeExpansion'
        val=fluid.enthalpy(pa,Ta)-pa./fluid.density(pa,Ta);
        u=-full(diag(val.jac{1}))./full(diag(val.jac{2}));
    case 'adiabaticExpansion'
        %val=fluid.enthalpy(pa,Ta)-pa./double(fluid.density(pa,Ta));
        val=fluid.enthalpy(pa,Ta)-pa./fluid.density(pa,Ta)  + double(pa).*(1./fluid.density(pa,Ta));
        u=-full(diag(val.jac{1}))./full(diag(val.jac{2}));
    case 'c_p'
        val=fluid.enthalpy(pa,Ta);
        u=full(val.jac{2});
    case 'comp_p'
        val=fluid.density(pa,Ta);
        u=full(val.jac{1})/double(val);
    case 'comp_t'
        val=fluid.density(pa,Ta);
        u=-full(val.jac{2})./double(val);
    otherwise
        error('no such property implemented')
end
end