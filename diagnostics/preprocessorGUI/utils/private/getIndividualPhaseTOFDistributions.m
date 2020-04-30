function [sphase,tof] = getIndividualPhaseTOFDistributions(data,tof,tof_ix,inj_ix,D,phase)

        itr = D.itracer(tof_ix, :);
        itr = [itr 1-sum(itr,2)];
        sphase = bsxfun(@times, data(:, phase), itr(:, inj_ix));
        sphase = cumsum(sphase);

end