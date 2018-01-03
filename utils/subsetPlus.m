function u = subsetPlus(u, v, subs)
    u(subs, :) = u(subs, :) + v;
end