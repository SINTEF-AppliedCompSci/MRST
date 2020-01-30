function pts = ellipticSurface3D(center, major_axis, minor_axis, major_axis_angle, ...
                 strike_angle, dip_angle, varargin)
    %{
   Create an elliptic surface that is apparoximated by a polygon. Addapted 
    from https://github.com/pmgbergen/porepy

    The rotation of the plane is calculated using three angles. First, the
    rotation of the major axis from the x-axis. Next, the surface is
    inclined by specifying the strike angle (which gives the rotation
    axis) measured from the x-axis, and the dip angle. All angles are
    measured in radians.

    Parameters:
    center (np.ndarray, size 3x1): Center coordinates of surface.
    major_axis (double): Length of major axis (radius-like, not
    diameter).
    minor_axis (double): Length of minor axis. There are no checks on
    whether the minor axis is less or equal the major.
    major_axis_angle (double, radians): Rotation of the major axis from
    the x-axis. Measured before strike-dip rotation, see above.
    strike_angle (double, radians): Line of rotation for the dip.
    Given as angle from the x-direction.
    dip_angle (double, radians): Dip angle, i.e. rotation around the
    strike direction.
    num_points (int, optional): Number of points used to approximate
    the ellipsis. Defaults to 16.

    Example:
    Surface centered at [0, 1, 0], with a ratio of lengths of 2,
    rotation in xy-plane of 45 degrees, and an incline of 30 degrees
    rotated around the x-axis.
    >>> frac = EllipticFracture(np.array([0, 1, 0]), 10, 5, np.pi/4, 0,
                            np.pi/6)

    %}
    opt = struct('numPoints', 16);
   
    opt = merge_options(opt, varargin{:});
    num_points = opt.numPoints;

    % First, populate polygon in the xy-plane
    angs = linspace(0, 2 * pi, num_points + 1);
    angs = angs(1:end-1);
    x = major_axis * cos(angs);
    y = minor_axis * sin(angs);
    z = zeros(1, numel(angs));
    ref_pts = [x', y', z'];

    % Rotate reference points so that the major axis has the right
    % orientation
    major_axis_rot = rotationMatrixFromVector(major_axis_angle, [0, 0, 1]);
    rot_ref_pts = ref_pts * major_axis_rot';

    % Then the dip
    % Rotation matrix of the strike angle
    strike_rot = rotationMatrixFromVector(strike_angle, [0, 0, 1]);
    % Compute strike direction
    strike_dir = [1, 0, 0] * strike_rot';
    dip_rot = rotationMatrixFromVector(dip_angle, strike_dir);

    dip_pts = rot_ref_pts* dip_rot';

    % Set the points
    pts = center + dip_pts;
