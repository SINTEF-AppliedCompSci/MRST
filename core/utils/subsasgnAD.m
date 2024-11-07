function u = subsasgnAD(u, ind, v)
% make sure that input u get converted to AD when v is AD
    if ~isa(u, 'ADI') && isa(v,'ADI') % u is a vector
        u = double2ADI(u, v);
    end
    u(ind) = v;
end