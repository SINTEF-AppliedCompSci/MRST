function writeTestsTAP_YAMLISH(tests, filename)
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