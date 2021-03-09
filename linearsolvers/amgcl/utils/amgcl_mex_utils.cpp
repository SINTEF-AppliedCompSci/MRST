#define GET_STRUCT_SCALAR(pa, fld) (mxGetScalar(mxGetField(pa, 0, fld)))

/* Relaxation */
struct relax_opts {
    int relax_id;
    int ilut_p;
    double ilut_tau;
    int iluk_k;
    double ilu_damping;
    double jacobi_damping;
    int chebyshev_degree;
    double chebyshev_lower;
    int chebyshev_power_iters;
};

const char * combine(std::string prefix, std::string suffix){
    std::string out = prefix + suffix;
    return out.c_str();
}

void setRelaxationStructMex(relax_opts &opt, const mxArray * pa, std::string prefix){
    /* Convert mex struct pointer to struct for amg */
    std::string tmp;
    tmp = prefix + "relaxation";
    opt.relax_id = (int)GET_STRUCT_SCALAR(pa,tmp.c_str());
    tmp = prefix + "ilut_p";
    opt.ilut_p = (int)GET_STRUCT_SCALAR(pa,tmp.c_str());
    tmp = prefix + "ilut_tau";
    opt.ilut_tau = GET_STRUCT_SCALAR(pa,tmp.c_str());
    tmp = prefix + "iluk_k";
    opt.iluk_k = (int)GET_STRUCT_SCALAR(pa,tmp.c_str());
    tmp = prefix + "ilu_damping";
    opt.ilu_damping = GET_STRUCT_SCALAR(pa,tmp.c_str());
    tmp = prefix + "jacobi_damping";
    opt.jacobi_damping = GET_STRUCT_SCALAR(pa,tmp.c_str());
    tmp = prefix + "chebyshev_degree";
    opt.chebyshev_degree = (int)GET_STRUCT_SCALAR(pa,tmp.c_str());
    tmp = prefix + "chebyshev_lower";
    opt.chebyshev_lower = GET_STRUCT_SCALAR(pa,tmp.c_str());
    tmp = prefix + "chebyshev_power_iters";
    opt.chebyshev_power_iters = (int)GET_STRUCT_SCALAR(pa,tmp.c_str());
}


void setRelaxationAMGCL(boost::property_tree::ptree & prm, std::string prefix, relax_opts opts){
    std::string relaxType = prefix + "type";
    switch(opts.relax_id) {
        case 1:
            prm.put(relaxType,  amgcl::runtime::relaxation::spai0);
            break;
        case 2:
            prm.put(relaxType,  amgcl::runtime::relaxation::gauss_seidel);
            break;
        case 3:
            prm.put(relaxType,  amgcl::runtime::relaxation::ilu0);
            prm.put(prefix + "damping", opts.ilu_damping);
            break;
        case 4:
            prm.put(relaxType,  amgcl::runtime::relaxation::iluk);
            prm.put(prefix + "damping", opts.ilu_damping);
            prm.put(prefix + "k", opts.iluk_k);
            break;
        case 5:
            prm.put(relaxType,  amgcl::runtime::relaxation::ilut);
            prm.put(prefix + "damping", opts.ilu_damping);
            prm.put(prefix + "p", opts.ilut_p);
            prm.put(prefix + "tau", opts.ilut_tau);
            break;
        case 6:
            prm.put(relaxType,  amgcl::runtime::relaxation::damped_jacobi);
            prm.put(prefix + "damping", opts.jacobi_damping);
            break;
        case 7:
            prm.put(relaxType,  amgcl::runtime::relaxation::spai1);
            break;
        case 8:
            prm.put(relaxType,  amgcl::runtime::relaxation::chebyshev);
            prm.put(prefix + "degree", opts.chebyshev_degree);
            prm.put(prefix + "lower", opts.chebyshev_lower);
            prm.put(prefix + "power_iters", opts.chebyshev_power_iters);
            break;
        default : mexErrMsgTxt("Unknown relax_id.");
    }
}

/* Coarsening */
struct amg_opts {
    int coarsen_id;
    int coarse_enough;
    bool direct_coarse;
    int max_levels;
    int ncycle;
    int npre;
    int npost;
    int pre_cycles;
    double aggr_eps_strong;
    double aggr_over_interp;
    double aggr_relax;
    double rs_eps_strong;
    double rs_eps_trunc;
    bool rs_trunc;
};

void setCoarseningStructMex(amg_opts &c_opt, const mxArray * pa){
    /* Convert mex struct pointer to struct for amg */
    c_opt.coarsen_id = (int)GET_STRUCT_SCALAR(pa,"coarsening");
    c_opt.coarse_enough = (int)GET_STRUCT_SCALAR(pa,"coarse_enough");
    c_opt.direct_coarse = GET_STRUCT_SCALAR(pa,"direct_coarse");
    c_opt.max_levels = (int)GET_STRUCT_SCALAR(pa,"max_levels");
    /* Define cycle */
    c_opt.ncycle = (int)GET_STRUCT_SCALAR(pa,"ncycle");
    c_opt.npre = (int)GET_STRUCT_SCALAR(pa,"npre");
    c_opt.npost = (int)GET_STRUCT_SCALAR(pa,"npost");
    c_opt.pre_cycles = (int)GET_STRUCT_SCALAR(pa,"pre_cycles");
    /* Coarsening options for general aggregation */
    c_opt.aggr_eps_strong = GET_STRUCT_SCALAR(pa,"aggr_eps_strong");
    /* Regular aggregation */
    c_opt.aggr_over_interp = GET_STRUCT_SCALAR(pa,"aggr_over_interp");
    /* Smoothed aggregation */
    c_opt.aggr_relax = GET_STRUCT_SCALAR(pa,"aggr_relax");
    /* Coarsening options for Ruge-Stuben coarsening */
    c_opt.rs_eps_strong = GET_STRUCT_SCALAR(pa,"rs_eps_strong");
    c_opt.rs_trunc = GET_STRUCT_SCALAR(pa,"rs_trunc");
    c_opt.rs_eps_trunc = GET_STRUCT_SCALAR(pa,"rs_eps_trunc");
}

void setCoarseningAMGCL(boost::property_tree::ptree & prm, std::string prefix, amg_opts options){
    std::string coarsetype = prefix + "coarsening.type";
    switch(options.coarsen_id) {
        case 1:
            prm.put(coarsetype,  amgcl::runtime::coarsening::smoothed_aggregation);
            prm.put(prefix + "coarsening.relax", options.aggr_relax);
            break;
        case 2:
            prm.put(coarsetype,  amgcl::runtime::coarsening::ruge_stuben);
            prm.put(prefix + "coarsening.eps_strong", options.rs_eps_strong);
            prm.put(prefix + "coarsening.do_trunc", options.rs_trunc);
            prm.put(prefix + "coarsening.eps_trunc", options.rs_eps_trunc);
            break;
        case 3:
            prm.put(coarsetype,  amgcl::runtime::coarsening::aggregation);
            prm.put(prefix + "coarsening.over_interp", options.aggr_over_interp);
            break;
        case 4:
            prm.put(coarsetype,  amgcl::runtime::coarsening::smoothed_aggr_emin);
            break;
        default : mexErrMsgTxt("Unknown coarsen_id: " + options.coarsen_id);
    }
    if(options.coarsen_id != 2){
        prm.put(prefix + "coarsening.aggr.eps_strong", options.aggr_eps_strong);
    }
    /* When is a level coarse enough */
    if (options.coarse_enough >= 0){
        prm.put(prefix + "coarse_enough", options.coarse_enough);
    }
    /* Use direct solver for coarse sys */
    prm.put(prefix + "direct_coarse", options.direct_coarse);
    /* Max levels */
    if (options.max_levels >= 0){
        prm.put(prefix + "max_levels", options.max_levels);
    }
    /* Number of cycles */
    if (options.ncycle >= 0){
        prm.put(prefix + "ncycle", options.ncycle);
    }
    /* Pre cycles */
    if (options.npre >= 0){
        prm.put(prefix + "npre", options.npre);
    }
    /* Post cycles */
    if (options.npost >= 0){
        prm.put(prefix + "npost", options.npost);
    }
    /* Pre cycles (precond) */
    if (options.pre_cycles >= 0){
        prm.put(prefix + "pre_cycles", options.pre_cycles);
    }
}
/* Krylov solver */
struct solver_opts {
    int solver_id;
    int L;
    int M;
    int K;
    int S;
    double delta;
    double omega;
    bool convex;
    bool always_reset;
    bool store_Av;
    bool replace;
};

void setSolverStructMex(solver_opts &opt, const mxArray * pa){
    opt.solver_id = (int)GET_STRUCT_SCALAR(pa,"solver");
    opt.L = (int)GET_STRUCT_SCALAR(pa,"bicgstabl_l");
    opt.M = (int)GET_STRUCT_SCALAR(pa,"gmres_m");
    opt.K = (int)GET_STRUCT_SCALAR(pa,"lgmres_k");
    opt.S = (int)GET_STRUCT_SCALAR(pa,"idrs_s");
    opt.delta = GET_STRUCT_SCALAR(pa,"bicgstabl_delta");
    opt.omega = GET_STRUCT_SCALAR(pa,"idrs_omega");
    opt.convex = GET_STRUCT_SCALAR(pa,"bicgstabl_convex");
    opt.always_reset = GET_STRUCT_SCALAR(pa,"lgmres_always_reset");
    opt.store_Av = GET_STRUCT_SCALAR(pa,"lgmres_store_av");
    opt.replace = GET_STRUCT_SCALAR(pa,"idrs_replacement");
}


void setSolverAMGCL(boost::property_tree::ptree & prm, std::string prefix, solver_opts options){
    std::string solvertype = prefix + "type";
    switch(options.solver_id) {
        case 1:
            prm.put(solvertype,  amgcl::runtime::solver::bicgstab);
            prm.put(prefix + "check_after", true);
            break;
        case 2:
            prm.put(solvertype,  amgcl::runtime::solver::cg);
            break;
        case 3:
            prm.put(solvertype,  amgcl::runtime::solver::bicgstabl);
            {
                prm.put(prefix + "L", options.L);
                prm.put(prefix + "delta", options.delta);
                prm.put(prefix + "convex", options.convex);
            }
            break;
        case 4:
            prm.put(solvertype,  amgcl::runtime::solver::gmres);
            {
                prm.put(prefix + "M", options.M);
            }
            break;
        case 5:
            prm.put(solvertype,  amgcl::runtime::solver::lgmres);
            {
                prm.put(prefix + "M", options.M);
                prm.put(prefix + "K", options.K);
                prm.put(prefix + "always_reset", options.always_reset);
                prm.put(prefix + "store_Av", options.store_Av);
            }
            break;
        case 6:
            prm.put(solvertype,  amgcl::runtime::solver::fgmres);
            {
                prm.put(prefix + "M", options.M);
            }
            break;
        case 7:
            prm.put(solvertype,  amgcl::runtime::solver::idrs);
            {
                prm.put(prefix + "s", options.S);
                prm.put(prefix + "omega", options.omega);
                prm.put(prefix + "replacement", options.replace);
            }
            break;
        default : mexErrMsgTxt("Unknown solver_id.");
    }
}
