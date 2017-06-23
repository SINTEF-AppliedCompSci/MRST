function [fx, cx] = splitFaceCellValue(flag, x, N, sz, upstr)
    nf = sz(1);
    nc = sz(2);
    
    if isa(x, 'ADI')
        n = numval(x);
    else
        n = numel(x);
    end
    
    switch n
        case nc
            % Cell-wise values only, use upstream weighting
            fx = upstr(flag, x);
            cx = x;
        case nf + nc
            % Face values first, then cell values
            fx = x(1:nf);
            cx = x((nf+1):end);
        case 2*nf + nc
            % Half face values
            subs = (1:nf)' + ~flag.*nf;
            fx = x(subs);
            cx = x((2*nf+1):end);
        case nf
            fx = x(1:nf);
            error('Not implemented yet');
        case 2*nf
            error('Not implemented yet');
        otherwise
            error('Did not find expected dimension of input');
    end
end