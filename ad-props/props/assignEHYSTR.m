function f = assignEHYSTR(f, ehystr, reg)
%
% Assigns flags for hysteresis computation, as well as input options for
% hysteresis model according to ECLIPSE EHYSTR keyword. 
% 

% flags for hysteresis
if strcmp(ehystr{5}, 'KR') || strcmp(ehystr{5}, 'BOTH')
    f.krHyst = 1;
elseif strcmp(ehystr{5}, 'PC') || strcmp(ehystr{5}, 'BOTH')
    f.pcHyst = 1;
end

% input options
f.ehystr = ehystr;

end