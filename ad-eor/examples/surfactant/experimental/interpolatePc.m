function [s, pc] = interpolatePc()
% Values from paper from Evans, Numerical Simulation of Water flooding Pancake Oil Column (1970)

%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}

    pc = [[.15  	46.00 ];
          [.18  	31.00 ];
          [.21  	22.86 ];
          [.25  	14.29 ];
          [.28  	9.16  ];
          [.30  	6.54  ];
          [.33  	3.92  ];
          [.35  	2.79  ];
          [.36  	2.47  ];
          [.38  	2.10  ];
          [.385	2.03  ];
          [.39  	1.99  ];
          [.41  	1.91  ];
          [.45  	1.81  ];
          [.50  	1.70  ];
          [.55  	1.60  ];
          [.60  	1.51  ];
          [.65  	1.44  ];
          [.70  	1.38  ];
          [.75  	1.32  ];
          [.80  	1.25  ];
          [.81  	1.23  ];
          [.85  	1.15  ];
          [.90  	0.35  ];
          [.94  	-2.24 ];
          [.96  	-5.08 ];
          [.97  	-7.18 ];
          [.98  	-10.10];
          [.99  	-14.00];
          [.9999 -31.00]];
    pc(:, 2) = convertFrom(pc(:, 2), psia);
    % Input deck values for capillary pressure are in bars.
    pc(:, 2) = convertTo(pc(:, 2), barsa);
    s = (0:0.01:1)';
    pc = interp1(pc(:, 1), pc(:, 2), s, 'spline');
end
