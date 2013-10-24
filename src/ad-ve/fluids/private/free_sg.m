function fsg = free_sg(sg,sGmax,opt)
        ineb=sg>=sGmax;
        %sg_free=sg-(loc_opt.sGmax-sg)*opt.res_gas;                
        fsg=((1-opt.res_oil)*sg-(sGmax*opt.res_gas))./(1-opt.res_gas-opt.res_oil);
        %fsg=(sg-(sGmax*opt.res_gas))./(1-opt.res_gas);
        %fsg(sg>=sGmax)=sg(sg>=sGmax);
        %assert(all(fsg>=-sqrt(eps)));
        fsg(ineb)=sg(ineb);
        fsg(fsg<0)=0.0*fsg(fsg<0);
end 