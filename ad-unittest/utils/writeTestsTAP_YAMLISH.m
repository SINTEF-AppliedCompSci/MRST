function writeTestsTAP_YAMLISH(tests, filename)
%Undocumented Utility Function

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

    nt = numel(tests);
    h = fopen(filename, 'w');
    if h == -1
        error('Unable to open file');
    end
    fprintf(h, '1..%d\n', nt);
    
    for i = 1:nt
        test = tests(i);
        if test.Passed
            okstr = 'ok';
        else
            okstr = 'not ok';
        end
        n = strsplit(test.Name, '/');
        
        fprintf(h, '%s %d - %s', okstr, i, n{1});
        if numel(n) > 1
            fprintf(h, '/%s', n{2});
        end
        fprintf(h, '\n');
        fprintf(h, '  ---\n');
        fprintf(h, '  duration_ms: %f\n', test.Duration);
        fprintf(h, '  ...\n');
    end
    fclose(h);
end
