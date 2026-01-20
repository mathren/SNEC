import numpy as np
import matplotlib as mpl
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
from plot_aux import SNEC_output_parser
import astropy.units as u


def plot_vel_radius_at_time_t(t, vel_out, ax=None, fig_name=None, **kwargs):
    data = SNEC_output_parser(vel_out)
    keys = np.array(list(data.keys()))
    times = keys * u.s
    try:
        units = t.unit
    except AttributeError:
        t *= u.s
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



if __name__ == "__main__":
    root = "/home/mrenzo/Documents/Research/codes/SNEC-1.01/Data/baseline/"
    answer = input("Want to plot the files in "+root+" ?[Y/n]")
    if answer.lower() != "y":
        root = input("Input folder?")
        print("Plotting files in "+root)
    vfile = root+"vel.xg"
    data = SNEC_output_parser(vfile)
    keys = np.array(list(data.keys()))
    times = keys * u.s
    colors = plt.cm.viridis(np.linspace(0,1, len(times)))
    fig = plt.figure()
    gs = gridspec.GridSpec(150, 100)
    ax = fig.add_subplot(gs[:, :])
    for i, t in enumerate(times):
        m, v = plot_vel_radius_at_time_t(t, vfile, ax=ax,
                                         c=colors[i], label=f"{t.to(u.h):.1f}")
    ax.set_xlabel(r"Initial mass coordinate $[M_{\odot}]$")
    ax.set_ylabel(r"$v \ [\mathrm{km\ s^{-1}}]$")
    plt.show()
