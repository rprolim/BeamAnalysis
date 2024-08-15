def healpix_to_2d(data, width, height, dec_range, ra_range):
    """
    Returns a 2D projection of an FITS file in (num_bands, height, width).

    Parameters:
    - data: FITS file, HEALPix map
        Desired data already opened.
    - width: int
        Desired width in pixels of projected map.
    - height: int
        Desired height in pixels of projected map.
    - dec_range: array-like
        Range of declination in degrees where dec_range[0] is the minimum and dec_range[1] is the maximum.
    - ra_range: array-like
        Range of right ascension in degrees where ra_range[0] is the minimum and ra_range[1] is the maximum.

    Returns:
    - projection: array-like
        2D projection of the FITS file in (num_bands, height, width).
    """
    if data.ndim > 1:
        nside = hp.get_nside(data[0])
        num_bands = data.shape[0]
    else:
        nside = hp.get_nside(data)
        num_bands = 1

    # Convert RA and Dec to theta (colatitude) and phi (longitude) in radians
    theta_range = np.radians([90.0 - dec_range[1], 90.0 - dec_range[0]])  # theta_min, theta_max
    phi_range = np.radians(ra_range)  # phi_min, phi_max

    # Create an array for the 2D projection
    projection = np.zeros((num_bands, height, width))

    # Set up the Mollweide projection grid
    theta = np.linspace(theta_range[0], theta_range[1], height)  # Latitude from 0 to 180 degrees
    phi = np.linspace(phi_range[0], phi_range[1], width)  # Longitude from -180 to 180 degrees
    phi, theta = np.meshgrid(phi, theta)  # Create a mesh grid for theta and phi

    # Convert theta, phi to HEALPix pixel indices
    pix_indices = hp.ang2pix(nside, theta, phi)

    if data.ndim > 1:
        for band in range(num_bands):
            projection[band, :, :] = data[band][pix_indices]  # Map HEALPix data to 2D grid
    else:
        projection[:, :] = data[pix_indices]

    return projection

def smooth_map(healpix_map, nbands=30, nu_min=980, nu_max=1260, D=40, in_degree=True, type_='gaussian'):
    """
    Process a HEALPix map by adjusting frequency bands, calculating beam window function,
    and smoothing the map.
    The values of this function are adjusted for the BINGO radio-telescope.

    Parameters:
    - healpix_map: array-like
        The HEALPix map data.
    - nbands: int
        Number of frequency bands.
    - nu_min: float
        Minimum frequency of observation.
    - nu_max: float
        Maximum frequency of observation.
    - D: float
        Diameter for the dish.
    - in_degree: bool
        Whether to use degrees in the beam window function.
    - type_: str
        The type of the beam window function.

    Returns:
    - g_smoothed_map: array-like
        The smoothed HEALPix map.
    """
    # Adjusting nu
    nu_step = (nu_max - nu_min) / nbands
    nu = np.around(np.arange(nu_min, nu_max, nu_step), decimals=2)

    # Adjusting ell
    nside = hp.get_nside(healpix_map)
    lmax = 3 * nside
    l = np.arange(lmax + 1)

    # FWHM (Delta theta_r)
    fwhm = model.fwhm_modelling(nu=nu[0], type_='smooth', D=D, in_degree=in_degree)

    # Adjusting theta
    theta = np.arange(0, 10, 0.01)

    # Get beam window function
    bl = model.bl_function(type_=type_, fwhm=fwhm, lmax=lmax, theta_=theta, input_unit='degree')

    # Visualize results
    smoothed_map = hp.smoothing(healpix_map, beam_window=bl)

    return smoothed_map

def array_to_healpix(projection, dec_range, ra_range, nside):
    """
    Convert a 3D array projection into a multi-band HEALPix map.

    Parameters:
    - projection: array-like 
        3D array of shape (num_bands, height, width).
    - dec_range: array-like
        Range of declination in degrees where dec_range[0] is the minimum and dec_range[1] is the maximum.
    - ra_range: array-like
        Range of right ascension in degrees where ra_range[0] is the minimum and ra_range[1] is the maximum.
    - nside: int
        HEALPix nside parameter.

    Returns:
    healpix_map: HEALPix map
        HEALPix map array of shape (num_bands, 12 * nside**2).
    """
    num_bands, height, width = projection.shape
    
    # Initialize HEALPix map array
    if num_bands > 1:
        healpix_map = np.full((num_bands, 12 * nside**2), hp.UNSEEN)
    else:
        healpix_map = np.full(12 * nside**2, hp.UNSEEN)

    # Convert RA and Dec to theta (colatitude) and phi (longitude) in radians
    theta_range = np.radians([90.0 - dec_range[1], 90.0 - dec_range[0]])  # theta_min, theta_max
    phi_range = np.radians(ra_range)  # phi_min, phi_max

    # Set up the Mollweide projection grid
    theta = np.linspace(theta_range[0], theta_range[1], height)  # Latitude from dec_min to dec_max
    phi = np.linspace(phi_range[0], phi_range[1], width)  # Longitude from ra_min to ra_max
    phi, theta = np.meshgrid(phi, theta)  # Create a mesh grid for theta and phi

    # Convert theta, phi to HEALPix pixel indices
    pix_indices = hp.ang2pix(nside, theta.flatten(), phi.flatten())

    if num_bands > 1:
        # Map values to HEALPix pixels for each band
        for band in range(num_bands):
            healpix_map[band, pix_indices] = projection[band].flatten()
    else:
        healpix_map[pix_indices] = projection.flatten()

    return healpix_map