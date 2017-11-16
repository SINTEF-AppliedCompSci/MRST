function dp = nozzleValve(v, rho, D, dischargeCoeff, flowtype)
% Nozzle valve model
A = pi*(D/2).^2;
if nargin < 5
    flowtype = 'massRate';
end
switch flowtype
    case 'velocity'
        % do nothing
    case 'volumeRate'
        v = v./A;
    case 'massRate'
        v = v./(rho*A);
    otherwise
        error(['Unknown flow type: ', flowtype]);
end

dp = - (sign(double(v)).*rho.*v.^2)./(2*dischargeCoeff.^2);
end
