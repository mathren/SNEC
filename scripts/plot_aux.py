import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.gridspec as gridspec
from MESAreader import get_src_col, Lsun, Rsun_cm, Msun, G_cgs, clight
import numpy as np
import re
import  astropy.units as u



class MultiColorPatch:
    def __init__(self, colors, alpha=1.0):
        self.colors = colors
        self.alpha = alpha


# Handler class that matplotlib will use
class MultiColorHandler:
    def legend_artist(self, legend, orig_handle, fontsize, handlebox):
        width = handlebox.width
        height = handlebox.height
        colors = orig_handle.colors
        alpha = orig_handle.alpha
        n_colors = len(colors)

        patches = []
        for i, color in enumerate(colors):
            patch = mpl.patches.Rectangle(
                [width * i / n_colors, 0],
                width / n_colors,
                height,
                facecolor=color,
                edgecolor='none',
                alpha=alpha,
                transform=handlebox.get_transform()
            )
            handlebox.add_artist(patch)
            patches.append(patch)

        return patches[0]


def HRD(fname, ax, annotate_radii=None, interp=None, scatter=False,  **kwargs):
    """
    make an HRD given the path fname to a history.data file on the provided ax
    annotate_radii = None or list of float values, if not None mark lines of constant radii for each provided value
    interp = None or float, if not None plot one point every `interp` interval of time
    kwargs = other plotting args
    """
    src, col = get_src_col(fname)
    logT = src[:, col.index("log_Teff")]
    logL = src[:, col.index("log_L")]
    if interp:
        t = src[:, col.index("star_age")]
        t_new = np.arange(0, max(t), interp)
        logL_interp = np.interp(t_new, t, logL)
        logT_interp = np.interp(t_new, t, logT)
        try:
            color = kwargs['c']
        except:
            try:
                color = kwargs['color']
            except:
                color='C0'
        ax.plot(logT, logL, lw=1, c=color, alpha=0.8)
        ax.scatter(logT_interp, logL_interp, **kwargs)
    else:  # not interpolating
        if not scatter:
            ax.plot(logT, logL, **kwargs)
        else:
            ax.scatter(logT, logL, **kwargs)
    if annotate_radii:
        annotate_radii_hrd(ax, radii=annotate_radii)
    return logT, logL


def annotate_radii_hrd(ax, radii=np.logspace(0, 3, base=10)):
    """
    give the axis object for an HRD plot (assumed to be in log10 Lsun and log10 Teff),
    and a list of radii in Rsun units, plots radii.
    Parameters:
    ----------
    ax: `mpl.ax` matplotlib axis object
    radii: `np.array`, optional, radii to mark
    """
    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()
    Tmax = 10.0 ** xmax
    Tmin = 10.0 ** xmin
    teff = np.linspace(Tmin, Tmax, 10)
    x = np.log10(teff)
    for r in radii:
        l = get_L_from_r_teff(r, teff)
        y = np.log10(l)
        ax.plot(x, y, c="#808080", ls=":", lw=1)
        # ax.text(x[5], y[5], f"{r:.0f}"+r"$\,R_\odot$", fontsize=20,
        # transform=ax.transData, zorder=0, c="#808080"),
        # rotation=np.((max(y)-min(y))/(max(x)-min(x))))
    # reset ylim
    ax.set_ylim(ymin, ymax)


def get_L_from_r_teff(radius, teff):
    # to annotate radii on HRD
    # Stephan Boltzman constant
    boltzm = 1.380649e-16  # cgs
    hbar = 6.62607015e-27 / (2 * np.pi)
    clight = 2.99792458e10
    sigma = (np.pi * np.pi * boltzm * boltzm * boltzm * boltzm) / (
        60 * hbar * hbar * hbar * clight * clight
    )
    # convert r to cm
    radius *= Rsun_cm
    # assume teff is in K
    l = 4.0 * np.pi * radius * radius * sigma * teff ** 4.0
    # convert to Lsun
    l = l / Lsun
    return l


def plot_Tc_rhoc(hfile, ax, **kwargs):
    src, col = get_src_col(hfile)
    logTc = src[:, col.index("log_center_T")]
    logRhoc = src[:, col.index("log_center_Rho")]
    ax.plot(logRhoc, logTc, **kwargs)
    return logTc, logRhoc


def get_iRLOF(bin_hfile):
    src, col = get_src_col(bin_hfile)
    rl_relative_overflow_1 = src[:, col.index("rl_relative_overflow_1")]
    iRLOF = rl_relative_overflow_1 >0
    return iRLOF


def calculate_lso_radius_and_spin_parameter(mass_bh, spin_parameter):
    """
    Calculate the radius and specific angular momentum for the marginally stable circular orbit
    around a Kerr black hole.

    Args:
        mass_bh (float): Mass of the black hole (in Msun).
        spin_parameter (float): Angular momentum a per unit mass (0 <= a <= 1).

    Returns:
        tuple: A tuple containing (r_lso_km, j_lso_cm2s).
            r_lso_cm (float): Radius of the marginally stable circular orbit (in km).
            j_lso_cm2s (float): Specific angular momentum (in cm^2/s).
    """
    # if not (0 <= spin_parameter <= 1):
    #     raise ValueError("Kerr parameter 'a' must be in the range [0, 1].")

    spin_parameter = abs(np.minimum(spin_parameter, 1.0))
    a = spin_parameter

    XX = (1.0 + a) ** (1.0 / 3.0) + (1.0 - a) ** (1.0 / 3.0)
    Z1 = 1.0 + XX * (1.0 - a ** 2.0) ** (1.0 / 3.0)
    Z2 = (3.0 * a ** 2.0 + Z1 ** 2.0) ** 0.5
    YY = ((3.0 - Z1) * (3.0 + Z1 + 2 * Z2)) ** 0.5

    r_lso1 = mass_bh * (3.0 + Z2 - YY)  # Prograde Orbit
    r_lso2 = mass_bh * (3.0 + Z2 + YY)  # Retrograde Orbit

    r_lso = r_lso1  # Only consider prograde (direct) orbits

    j_lso = r_lso * r_lso * mass_bh ** 0.5 / (r_lso ** (3.0 / 2.0) + a * mass_bh ** (3.0 / 2.0))
    j_lso = j_lso * Msun * G_cgs / clight  # Convert to cm^2/s
    r_lso = r_lso1 * Msun * G_cgs / (clight * clight)  # Convert to cm

    return r_lso, j_lso


def get_BE_from_pfile(pfile):
    """Calculates the binding energy profile of the star. See Eq. 6 in
    Dewi & Tauris 2000 (but change sign, binding energy>0 if layer is bound).

    The binding energy is calculated integrating the potential.

    Parameters
    ----------
    pfile    : `MESA profile*.data` file
               assumed to contain the columns mass, energy, radius, and dm
    Returns
    -------
    BE       : `np.array` binding energy in cgs units
    """
    # from time import time
    # get the data in cgs units
    src, col = get_src_col(pfile)
    m = src[:, col.index("mass")] * Msun  # g
    try:
        dm = src[:, col.index("dm")]  # g
    except ValueError:
        dm = np.diff(m, append=0)
    r = src[:, col.index("radius")] * Rsun_cm  # cm
    psi = -1.0 * G_cgs * np.divide(m, r)  # erg
    # change sign: BE is the energy (>0) to provide to unbind the star
    BE = -1.0 * np.cumsum(np.multiply(psi, dm))
    return np.asarray(BE, dtype=float)


def find_CO_core_from_pfile(pfile,
                            co_core_boundary_he4_fraction=0.1,
                            min_boundary_fraction=0.1):
    """ define CO core as MESA:
    outermost location where he4 mass fraction is <= co_core_boundary_he4_fraction,
    and c12+o16 mass fraction >= min_boundary_fraction.

    Returns CO core mass in Msun units
    """
    src, col = get_src_col(pfile)
    c12 = src[:, col.index("c12")]
    o16 = src[:, col.index("o16")]
    he4 = src[:, col.index("he4")]
    m = src[:, col.index("mass")]
    co = c12+o16
    ihe_poor = (he4 <= co_core_boundary_he4_fraction)
    ico_rich = (co >= min_boundary_fraction)
    return max(m[ico_rich & ihe_poor])


def find_He_core_from_pfile(pfile,
                            he_core_boundary_h_fraction=0.1,
                            min_boundary_fraction=0.1):
    """ define He core as MESA:
    outermost location where h1 mass fraction is <= he_core_boundary_h_fraction,
    and he4 mass fraction >= min_boundary_fraction.

    Returns CO core mass in Msun units
    """
    src, col = get_src_col(pfile)
    h1 = src[:, col.index("h1")]
    he4 = src[:, col.index("he4")]
    m = src[:, col.index("mass")]
    ih_poor = (h1 <= he_core_boundary_h_fraction)
    ihe_rich = (he4 >= min_boundary_fraction)
    return max(m[ihe_rich & ih_poor])



def SNEC_output_parser(outfile):
    """
    Parse the file and return a dictionary with time values as keys
    and numpy arrays with two columns as values.

    Parameters:
    -----------
    filename : str
        Path to the input file

    Returns:
    --------
    data: np.array [time, mass, radius]
    """
    data_dict = {}

    with open(outfile, "r") as f:
        content = f.read()

    # Split by "Time =" to get each time block
    blocks = re.split(r'"Time =', content)

    for block in blocks[1:]:  # Skip the first empty split
        lines = block.strip().split('\n')

        # Extract time value from first line
        time_line = lines[0]
        time_value = float(time_line.split()[0])

        # Parse data lines
        data_lines = []
        for line in lines[1:]:
            line = line.strip()
            if line and (not line.startswith('"')) and (not line.startswith('#')):  # Skip empty lines and closing quotes
                parts = line.split()
                if len(parts) >= 2:
                    try:
                        col1 = float(parts[0].replace('E+0', 'E+'))
                        col2 = float(parts[1].replace('E+0', 'E+'))
                        data_lines.append([col1, col2])
                    except ValueError:
                        continue

        # Convert to numpy array
        if data_lines:
            data_dict[time_value] = np.array(data_lines)

    return data_dict




def plot_vel_radius_at_time_t(t, vel_out, ax=None, fig_name=None, **kwargs):
    data = SNEC_output_parser(vel_out)
    keys = np.array(list(data.keys()))
    times = keys * u.s
    try:
        units = times.unit
    except AttributeError:
        times *= u.s
    index_time_of_interest = np.argmin(np.absolute(times-t))
    key_of_interest = keys[index_time_of_interest]
    mass = data[key_of_interest][:, 0] * u.g
    vel = data[key_of_interest][:,1] * u.cm/u.s
    if not ax:
        fig = plt.figure()
        gs = gridspec.GridSpec(150, 100)
        ax = fig.add_subplot(gs[:, :])
    ax.plot(mass.to(u.Msun), vel.to(u.km/u.s), **kwargs)
    if fig_name:
        ax.set_ylabel(r"$v \ [\mathrm{km\ s^{-1}}]$")
        ax.set_xlabel(r"$M \ [M_{\odot}]$")
        plt.savefig(fig_name)
    return mass, vel



def plot_mass_radius(t, mass_out, ax=None, **kwargs):
    data = SNEC_output_parser(mass_out)
    keys = np.array(list(data.keys()))
    times = keys * u.s
    try:
        units = times.unit
    except AttributeError:
        times *= u.s
    index_time_of_interest = np.argmin(np.absolute(times-t))
    key_of_interest = keys[index_time_of_interest]
    mass = data[key_of_interest][:, 1] * u.g
    radius = data[key_of_interest][:, 0] * u.cm
    if not ax:
        fig = plt.figure()
        gs = gridspec.GridSpec(150, 100)
        ax = fig.add_subplot(gs[:, :])
    p = ax.scatter(mass.to(u.Msun), radius.to(u.cm), c=[t.to(u.h).value]*len(mass), **kwargs)
    return (mass, radius, p)



def plot_v_radius_time(t, vel_out, mass_out, ax=None,
                       annotate_inner_boundary=True,
                       annotate_pre_expl_R=True,
                       **kwargs):
    vel_data = SNEC_output_parser(vel_out)
    keys = np.array(list(vel_data.keys()))
    vel_times = keys * u.s
    try:
        units = t.unit
    except AttributeError:
        vel_times *= u.s
    mass_data = SNEC_output_parser(mass_out)
    keys = np.array(list(mass_data.keys()))
    mass_times = keys * u.s
    try:
        units = mass_times.unit
    except AttributeError:
        mass_times *= u.s
    # sanity check
    # print(vel_times == mass_times)
    index_time_of_interest = np.argmin(np.absolute(vel_times-t))
    key_of_interest = keys[index_time_of_interest]
    mass = mass_data[key_of_interest][:, 1] * u.g
    radius = mass_data[key_of_interest][:, 0] * u.cm
    vel = vel_data[key_of_interest][:, 1] * u.cm/u.s
    vel = vel.to(u.km/u.s)
    if not ax:
        fig = plt.figure()
        gs = gridspec.GridSpec(150, 100)
        ax = fig.add_subplot(gs[:, :])
    if annotate_inner_boundary:
        i_min_m = np.argmin(mass)
        ax.axvline(np.log10(radius[i_min_m].value), 0, 1, zorder=0,
                   ls='--', lw=1, c='k')
    if annotate_pre_expl_R:
        ax.axvline(np.log10(max(radius.value)), 0, 1, zorder=0, ls='--', lw=1)
    ax.plot(np.log10(radius.value), vel, **kwargs)
    # ax.scatter(np.log10(radius.value), vel, **kwargs)


def get_times(data_file):
    data = SNEC_output_parser(data_file)
    keys = np.array(list(data.keys()))
    times = keys * u.s
    return times
