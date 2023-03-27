#ifndef OPTIONS_H
#define OPTIONS_H

/**
 *  A simple options structure. Same structure as used by matlab.
 */
typedef struct
{
    /// Flag to say if we should be verbose and print everything out.
    bool    verbose;
    /// Flag to say if we should compute time step satisfying CFL.
    bool    computedt;
    /// Flag to say if we should integrate vertically.
    bool    intVert;
    /// Flag to say if we should integrate the porosity vertically.
    bool    intVert_poro;
    /// Should we solve with a semi-implicit solver?
    bool    semi_implicit;
    /// Check if we should use upwind gravity.
    bool    grav_upwind;
    /// Check if we should use dif.
    bool    no_dif;
    /// Check if it's a central scheme.
    bool    central;
    /// The current time step.
    double  dt;
    /// The tolerance for overflowing heights.
    double  heightWarn;
    /// Norm of gravity component.
    double  gravity;
    /// The flag
    int     flag;
    /// The well structure.
    Well<double>    *wells;
    /// The boundary conditions.
    BC<double>      *bc;

} opt_t;

#endif
