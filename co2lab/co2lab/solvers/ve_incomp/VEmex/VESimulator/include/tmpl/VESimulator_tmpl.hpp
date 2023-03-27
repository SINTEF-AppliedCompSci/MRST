template<class Real>
VESimulator<Real>::VESimulator()
{
    m_state         = NULL;
    m_grid          = NULL;
    m_rock          = NULL;
    m_fluid         = NULL;
    m_wells         = NULL;
    m_bc            = NULL;

    verbose         = false;
    computeDt       = false;
    intVert         = true;
    intVert_poro    = false;
    semi_implicit   = false;
    no_dif          = false;
    central         = false;
    gravity_upwind  = false;
    heightWarn      = (Real) 1.0e-16;
    gravity         = (Real) 0.0;
}

template<class Real>
VESimulator<Real>::VESimulator(State<Real> *state, Grid2D<Real> *grid,
        Rock<Real> *rock, Fluid<Real> *fluid, opt_t *opts)
{
    m_state         = state;
    m_grid          = grid;
    m_rock          = rock;
    m_fluid         = fluid;

    if(opts->wells) {
        m_wells         = new Well<Real>(opts->wells);
        delete opts->wells;
    } else {
        m_wells         = new Well<Real>();
    }

    if(opts->bc) {
        m_bc            = new BC<Real>(opts->bc);
        delete opts->bc;
    } else {
        m_bc            = new BC<Real>();
    }

    verbose         = opts->verbose;
    computeDt       = opts->computedt;
    intVert         = opts->intVert;
    intVert_poro    = opts->intVert_poro;
    semi_implicit   = opts->semi_implicit;
    no_dif          = opts->no_dif;
    central         = opts->central;
    gravity_upwind  = opts->grav_upwind;
    heightWarn      = (Real) opts->heightWarn;
    timestep        = (Real) opts->dt;
    gravity         = (Real) opts->gravity;
}

template<class Real>
VESimulator<Real>::~VESimulator()
{
    if(m_state)
        delete m_state;

    if(m_grid)
        delete m_grid;

    if(m_rock)
        delete m_rock;

    if(m_fluid)
        delete m_fluid;

    if(m_wells)
        delete m_wells;

    if(m_bc)
        delete m_bc;
}

template<class Real>
State<Real>* VESimulator<Real>::getState() const
{
    return m_state;
}

template<class Real>
void VESimulator<Real>::changeState(State<Real> *state)
{
    if(m_state)
        delete m_state;
    m_state =  state;
}

template<class Real>
void VESimulator<Real>::changeOptions(opt_t *opts)
{
    if(m_wells)
        delete m_wells;
    if(m_bc)
        delete m_bc;

    if(opts->wells) {
        m_wells = new Well<Real>(opts->wells);
        delete opts->wells;
    } else {
        m_wells = new Well<Real>();
    }

    if(opts->bc) {
        m_bc            = new BC<Real>(opts->bc);
        delete opts->bc;
    } else {
        m_bc            = new BC<Real>();
    }

    verbose         = opts->verbose;
    computeDt       = opts->computedt;
    intVert         = opts->intVert;
    intVert_poro    = opts->intVert_poro;
    semi_implicit   = opts->semi_implicit;
    no_dif          = opts->no_dif;
    central         = opts->central;
    gravity_upwind  = opts->grav_upwind;
    heightWarn      = (Real) opts->heightWarn;
    timestep        = (Real) opts->dt;
    gravity         = (Real) opts->gravity;
}

template<class Real>
Real* VESimulator<Real>::sources()
{
    Real *q = new Real[m_grid->cells->num];
    std::fill(q, q + m_grid->cells->num, 0.0);

    m_wells->contribute_source(q, m_grid->cells->H);

    return q;
}

// Transport loop for the CPU, both single and multithreaded.
template<class Real>
void VESimulator<Real>::doTimeStep(Real& tf)
{
    const int numthreads = omp_get_max_threads();
    // Set the number of threads to use
    Real dt = std::min(timestep, tf);
    assert(dt > 0 && tf > 0);

    Real *h = m_state->height;
    Real *h_max = m_state->max_height;

    // Get sources.
    Real *q = sources();
    Real *d_flux = m_state->flux;

    // Compute the magnitude.
    Real rho_diff = m_fluid->rho[0] - m_fluid->rho[1];

    // Integrated permeability from 0 to inf.
    std::vector<Real> kr_H( m_grid->cells->num );
    std::vector<Real> dz( m_grid->cells->num );

    // Gravity for each face.
    std::vector<Real> grav( m_grid->faces->num );

    // Average porosity for each cell.
    std::vector<Real> poro( m_grid->cells->num );

    // z-component of normalized normal vector for each cell.
    std::vector<Real> n( m_grid->cells->num );

    Real *perm = m_rock->permeability;
    Real *norm = m_grid->cells->normals;

    // Mobility for each cell.
    std::vector<Real> mob( 2*m_grid->cells->num );
    std::vector<Real> dmob( 2*m_grid->cells->num );
    //const size_t numFaces = (m_grid->faces->num + numthreads - 1)/numthreads;
    //const size_t numCells = (m_grid->cells->num + numthreads - 1)/numthreads;




    std::vector<Real> dtkr_H( m_grid->cells->num );



//#pragma omp parallel
    {
        // Compute gravity vector.
#pragma omp parallel for
       for(size_t idx=0; idx < m_grid->faces->num; ++idx){
                int c1 = m_grid->faces->neighbors[2*idx];
                int c2 = m_grid->faces->neighbors[2*idx+1];

                //z_diff[idx] = m_grid->cells->z[c1] - m_grid->cells->z[c2];
                Real z_diff = m_grid->cells->z[c1] - m_grid->cells->z[c2];

                Real cdiff[2];
                Real e[2];
                cdiff[0] = m_grid->cells->centroids[2*c1] -
                    m_grid->cells->centroids[2*c2];
                cdiff[1] = m_grid->cells->centroids[2*c1+1] -
                    m_grid->cells->centroids[2*c2+1];
                e[0] = m_grid->nodes->coords[2*m_grid->faces->nodes[2*idx]] -
                    m_grid->nodes->coords[2*m_grid->faces->nodes[2*idx+1]];
                e[1] = m_grid->nodes->coords[2*m_grid->faces->nodes[2*idx]+1] -
                    m_grid->nodes->coords[2*m_grid->faces->nodes[2*idx+1]+1];

                Real c_ij = sqrt(z_diff*z_diff + cdiff[0]*cdiff[0] +
                        cdiff[1]*cdiff[1]);
                Real e_ij = sqrt(e[0]*e[0] + e[1]*e[1]);
                grav[idx] = gravity*e_ij/c_ij;
        }

        // Integrate permeability and average the porosity in z direction, also
        // calculate the normalized normal and pore volume
#pragma omp parallel for
        for(size_t idx = 0; idx < poro.size(); idx++) {
                // Integrate verticaly and average porosity.
                int from = m_grid->cells->columnPos[idx];
                int to = m_grid->cells->columnPos[idx+1];
                int num = to - from;
                int cell3D;

                kr_H[idx] = 0.0;
                dtkr_H[idx] = 0.0;

                while(from < to) {
                    cell3D = m_grid->columns->cells[from];

                    poro[idx] +=m_rock->porosity[cell3D]*m_grid->columns->dz[from];
                    kr_H[idx] += perm[m_grid->columns->cells[from]] *
                        m_grid->columns->dz[from];
                    dtkr_H[idx] += perm[m_grid->columns->cells[from]];

                    from++;
                }
                poro[idx] /= m_grid->cells->H[idx];
                dtkr_H[idx] /= num;

                // z value of the normalized normal vector.
                Real sum = sqrt(norm[3*idx]*norm[3*idx] +
                        norm[3*idx+1]*norm[3*idx+1] +
                        norm[3*idx+2]*norm[3*idx+2]);
                n[idx] = norm[3*idx+2]/sum;

        }
    }
    if(computeDt){
        dt=this->computeTimeStep(dt, d_flux, q, rho_diff, poro, grav, dtkr_H);
    }else{
        dt=tf;
    }

    // Run transport loop.
    Real t    = 0;
    Real *new_h = new Real[m_grid->cells->num];
    while(t < tf) {
        // Compute cell mobilities
#pragma omp parallel for
        for(int i = 0; i < m_grid->cells->num; i++) {
            Real kr;
            Real dkr;
            m_grid->cellRelperm(kr ,dkr , i, h[i], perm);
            //m_grid->cellMobility(i, h[i], perm);
            //Real dkr = m_grid->dcellMobility(i, h[i], perm);

            mob[2*i] = m_fluid->kwm[0]*kr/m_fluid->mu[0];
            mob[2*i+1] = m_fluid->kwm[1]*(kr_H[i] - kr)/m_fluid->mu[1];
            dmob[2*i] = m_fluid->kwm[0]*dkr/m_fluid->mu[0];
            dmob[2*i+1] = m_fluid->kwm[1]*dkr/m_fluid->mu[1];
        }

        // Loop over cells and it's faces
#pragma omp parallel
        {
        Real dt_loc=dt;
#pragma omp for nowait
        for(int i = 0; i < m_grid->cells->num; i++) {
            // Calculate pore volume
            Real pv = poro[i]*m_grid->cells->volumes[i];
            dz[i] = 0.0;
            for(int j = 0; j < 4; j++) {
                int face = m_grid->cells->faces[4*i+j];

                int c1 = m_grid->faces->neighbors[2*face];
                int c2 = m_grid->faces->neighbors[2*face+1];
                if(c1 != -1 && c2 != -1) {
                    Real z_diff = m_grid->cells->z[c1] - m_grid->cells->z[c2];
                    Real h_diff = n[c1]*h[c1] - n[c2]*h[c2];
                    Real g_flux = -grav[face]*(z_diff + h_diff)*rho_diff;

                    Real faceMob[2];
                    Real dfaceMob[2];
                    // Compute face mobility.
                    if( (g_flux*d_flux[face] >= 0) ) {
                        if(d_flux[face] > 0) {
                            faceMob[0] = mob[2*c1];
                            dfaceMob[0] = dmob[2*c1];
                        } else {
                            faceMob[0] = mob[2*c2];
                            dfaceMob[0] = dmob[2*c2];
                        }

                        if(d_flux[face] - faceMob[0]*g_flux > 0) {
                            faceMob[1] = mob[2*c1 + 1];
                            dfaceMob[1] = dmob[2*c1 + 1];
                        } else {
                            faceMob[1] = mob[2*c2 + 1];
                            dfaceMob[1] = dmob[2*c2 + 1];
                        }
                    } else {
                        if(d_flux[face] > 0) {
                            faceMob[1] = mob[2*c1+1];
                            dfaceMob[1] = dmob[2*c1+1];
                        } else {
                            faceMob[1] = mob[2*c2+1];
                            dfaceMob[1] = dmob[2*c2+1];
                        }

                        if(d_flux[face] + faceMob[1]*g_flux > 0) {
                            faceMob[0] = mob[2*c1];
                            dfaceMob[0] = dmob[2*c1];
                        } else {
                            faceMob[0] = mob[2*c2];
                            dfaceMob[0] = dmob[2*c2];
                        }
                    }

                    Real tmob = faceMob[0] + faceMob[1];
                    Real fw_face = faceMob[0]/tmob;

                    Real diff = fw_face*(d_flux[face] + faceMob[1]*g_flux);
                    {
                        Real v_co2 = diff;
                        Real v_water = (1-fw_face)*(d_flux[face] - faceMob[1]*g_flux);
                        Real ff=faceMob[1]*dfaceMob[0]*std::abs(v_co2)/(tmob*faceMob[0])+
                                faceMob[0]*dfaceMob[1]*std::abs(v_water)/(tmob*faceMob[1])+
                            std::abs(grav[face]*rho_diff*(faceMob[1]*fw_face));
                        assert(ff>=0);
                        Real dt_tmp=pv/ff;
                        dt_tmp *= 0.5*(1-(m_fluid->sr+m_fluid->sw));
                        dt_loc=std::min(dt_loc,dt_tmp);
                     //dt=std::min(dt_tmp,dt);//gives nan
                    }


                    diff *= (c1 == i ? 1 : -1);
                    dz[i] += diff;
                } else {
                    Real faceMob[2];
                    if(d_flux[face] > 0) {
                        if(c1!=-1){
                            faceMob[1] = mob[2*c1+1];
                            faceMob[0] = mob[2*c1];
                        }else{
                            faceMob[1]=1;
                            faceMob[0]=0;
                        }
                    } else {
                        if(c2!=-1){
                            faceMob[1] = mob[2*c2+1];
                            faceMob[0] = mob[2*c2];
                        }else{
                            faceMob[1]=1;
                            faceMob[0]=0;
                        }
                    }
                    Real sum = faceMob[0] + faceMob[1];
                    Real fw_face = faceMob[0]/sum;

                    Real diff = fw_face*d_flux[face];
                    diff *= (c1 == i ? 1 : -1);
                    dz[i] += diff;
                }
            }

            Real f_w = mob[2*i]/(mob[2*i] + mob[2*i+1]);
            dz[i] -= (std::max(q[i],(Real) 0.0) + std::min(q[i], (Real) 0.0)*f_w);
        }
        // this is moved here to avoid syncoisation problems between treads
#pragma omp critical
        {
            dt = std::min(dt_loc,dt);
        }

}
//should not be nessesary
//#pragma omp barrier

#pragma omp parallel for
        for(int i = 0; i < m_grid->cells->num; i++) {
            Real pv = poro[i]*m_grid->cells->volumes[i];
            // the false block is keept to reproduce old results if wanted.
            if(false){
                if(h_max[i] > h[i]) {
                    pv *= (1.0 - (m_fluid->sr + m_fluid->sw));
                } else {
                    pv *= (1.0 - m_fluid->sw);
                }

                // Compute new height
                new_h[i] = h[i] - dt*dz[i]/pv;

                new_h[i] = std::min(new_h[i], m_grid->cells->H[i]);
                new_h[i] = std::max(new_h[i], (Real) 0.0);

                h_max[i] = std::max(new_h[i], h_max[i]);
            }else{
                Real pv = poro[i]*m_grid->cells->volumes[i];
                Real   old_vol = pv*(h[i]*(1.0 - m_fluid->sw)+(h_max[i]-h[i])*m_fluid->sr);
                Real   new_vol = old_vol-dt*dz[i];
                Real   max_vol = pv*h_max[i]*(1.0 - m_fluid->sw);
                //if(new_vol<0){
                //    std::cout << i << new_vol << max_vol;
                //}
                if(new_vol>max_vol){
                    h_max[i]=new_vol/(pv*(1-m_fluid->sw));

                    h_max[i]= std::min(h_max[i], m_grid->cells->H[i]);
                    //	     	    std::cout << h_max[i];
                }
                new_h[i]=(new_vol-h_max[i]*pv*m_fluid->sr)/(pv*(1-(m_fluid->sr+m_fluid->sw)));

                if (new_h[i] > m_grid->cells->H[i]){
                    if(verbose){
                        fprintf(stderr, "h = %f > H in cell: %i \n ",new_h[i], i);
                    }
                    new_h[i] = std::min(new_h[i],  m_grid->cells->H[i]);
                }
                if (new_h[i] < 0 ){
                    if(verbose){
                        fprintf(stderr, "h = %f < 0 in cell: %i \n ", new_h[i], i);
                    }
                    new_h[i] = std::max(new_h[i], (Real) 0.0);
                }
            }
        }
        // Swap old height for new height.
        memcpy(h, new_h, m_grid->cells->num*sizeof(Real));
       if(verbose){
            fprintf(stderr, "Time step is %f  [day]\n",dt/(60*60*24));
        }
        t += dt;

        dt = std::min(dt, tf - t);
    }

    delete[] q;
    delete[] new_h;

    m_state->height = h;
    m_state->max_height = h_max;
}
template<class Real>
Real VESimulator<Real>::computeTimeStep(const Real& dT,
                                        const Real* d_flux,
                                        const Real* q,
					const Real& rho_diff,
                                        const std::vector<Real>& poro,
					const std::vector<Real>& grav,
                                        const std::vector<Real>& dtkr_H){
    // Some variables for computing the timestep.
    Real   maxheight = m_grid->cells->H[0];
#pragma omp parallel shared(maxheight)
    {
        Real maxheight_loc=maxheight;
#pragma omp for nowait
    for(size_t i = 0; i < m_grid->cells->num; i++) {
        maxheight_loc=std::max(maxheight_loc,m_grid->cells->H[i]);
    }
#pragma omp critical
        {
        maxheight=std::max(maxheight,maxheight_loc);
        }
    }


 Real dt=0;
 Real gravdt = dT;
 Real fluxdt = dT;
 Real  dfgrav = 0.0;
 Real  dfflux = 0.0;
#pragma omp parallel shared(gravdt,fluxdt,dfgrav,dfflux)
{

	  Real add = maxheight/1000.0;

          Real  dfgrav_loc = 0.0;
          Real dfflux_loc = 0.0;

        // Compute gravity vector.
	const int numthreads = omp_get_max_threads();
        const size_t threadnum = omp_get_thread_num();
	const size_t numFaces = (m_grid->faces->num + numthreads - 1)/numthreads;
	const size_t numCells = (m_grid->cells->num + numthreads - 1)/numthreads;
        //    const size_t dthCells = (1001 + numthreads - 1)/numthreads;
        // Compute time step.
#pragma omp for nowait
        for(size_t idx = 1; idx < 1001; idx++) {
                Real dth[2];
                dth[0] = (idx-1)*add;
                dth[1] = idx*add;

                Real mob_avg[2];
                mob_avg[0] = m_fluid->kwm[0]*dth[0]/m_fluid->mu[0];
                mob_avg[1] = m_fluid->kwm[1]*(maxheight-dth[0])/m_fluid->mu[1];
                Real sum = mob_avg[0] + mob_avg[1];

                Real oldgrav = mob_avg[0]*mob_avg[1]/sum;
                Real oldflux = mob_avg[0]/sum;

                mob_avg[0] = m_fluid->kwm[0]*dth[1]/m_fluid->mu[0];
                mob_avg[1] = m_fluid->kwm[1]*(maxheight-dth[1])/m_fluid->mu[1];
                sum = mob_avg[0] + mob_avg[1];

                Real newgrav = mob_avg[0]*mob_avg[1]/sum;
                Real newflux = mob_avg[0]/sum;

                Real dfg = fabs((newgrav - oldgrav)/add);
                Real dff = fabs((newflux - oldflux)/add);
                dfgrav_loc = std::max(dfgrav_loc, dfg);
                dfflux_loc = std::max(dfflux_loc, dff);

        }
#pragma omp critical
        {
           dfgrav = std::max(dfgrav, dfgrav_loc);
           dfflux = std::max(dfflux, dfflux_loc);
        }
#pragma omp barrier
        Real fluxdt_loc=fluxdt;
        Real gravdt_loc=gravdt;
#pragma omp for nowait
        for(size_t idx = 0; idx < m_grid->cells->num; idx++) {
                Real grav_in = 0.0f;
                Real grav_out = 0.0f;
                Real flux_in = 0.0f;
                Real flux_out = 0.0f;

                for(int j = 0; j < 4; j++) {
                    int facenum = m_grid->cells->faces[4*idx + j];

                    int c1 = m_grid->faces->neighbors[2*facenum];
                    int c2 = m_grid->faces->neighbors[2*facenum + 1];
                    if(c1 != -1 && c2 != -1) {
                        Real z_diff = m_grid->cells->z[c1] -
                            m_grid->cells->z[c2];
                        Real harmean = 2.0/(1.0/dtkr_H[c1] + 1.0/dtkr_H[c2]);

                        Real gravflux = grav[idx]*rho_diff*z_diff*harmean;

                        if(gravflux < 0) {
                            if(c1 == idx) {
                                grav_in += gravflux;
                            } else {
                                grav_out += gravflux;
                            }
                        } else {
                            if(c1 == idx) {
                                grav_out += gravflux;
                            } else {
                                grav_in += gravflux;
                            }
                        }

                        if(d_flux[facenum] < 0) {
                            if(c1 == idx) {
                                flux_in += d_flux[facenum];
                            } else {
                                flux_out += d_flux[facenum];
                            }
                        } else {
                            if(c1 == idx) {
                                flux_out += d_flux[facenum];
                            } else {
                                flux_in += d_flux[facenum];
                            }
                        }
                    }
                }

                Real max_grav = std::max(fabs(grav_in), fabs(grav_out));
                Real max_flux = std::max(
                fabs(flux_in) + std::max(q[idx], (Real) 0.0),
                fabs(flux_out) - std::min(q[idx], (Real) 0.0) );

                Real pv = poro[idx]*m_grid->cells->volumes[idx];
                if(max_grav != 0.0) {
                    Real step = pv/max_grav/dfgrav;
                    gravdt_loc = std::min(gravdt_loc, step);

                }

                if(max_flux != 0.0) {
                   Real step = pv/max_flux/dfflux;
                   fluxdt_loc = std::min(fluxdt_loc, step);
                }
        }
#pragma omp critical
       {
        fluxdt = std::min(fluxdt,fluxdt_loc);
        gravdt = std::min(gravdt,gravdt_loc);
       }

#pragma omp barrier

#pragma omp single
        {
            dt  = std::min(fluxdt, gravdt);
            dt *= 0.4*(1.0 - (m_fluid->sr + m_fluid->sw));
        }
}


return dt;
}
