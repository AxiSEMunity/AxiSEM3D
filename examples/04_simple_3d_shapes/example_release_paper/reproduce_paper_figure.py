import numpy as np
import obspy
from obspy.core.trace import Trace
from obspy.core.stream import Stream
import matplotlib.pyplot as plt
import sys
sys.path.append('../src/')
from clrs import Iridescent

def norm(x):
    return x/np.max(np.abs(x))

# Converts obspy trace into x and y for mpl plotting
def obspy_gen_mpl(tr):
    x = np.linspace(0, tr.stats.npts*tr.stats.delta,  tr.stats.npts)
    y = tr.data
    return x,y

# Input parameters and paths:
save_figs_as_pdf = True         # Either saves as PDF of plt.show() to screen
ddir = './datas/station_1'      # Path to data Station_1 directory from simulation
odir = './Figures/'             # Path to output plots
net = 'YY'                      # Network label
nstns = 15                      # Number of stations on each azimuth

norm_traces = True              # Normalise record section trace amplitudes

# Filter frequencies
fmin = 1                      # Minimum frequency in Hz
fmax = 7.5                      # Maximum frequency in Hz

# Load hex codes for a sequential, colour-blind-safe colour scheme
hex = Iridescent()


# Theoretical P wave arrival time:
arr_time = ((np.arange(1,16)*2000)**2 + 7500**2)**0.5 / 1500

# Load time data only once -- shared by all time-series
time = np.loadtxt(f"{ddir}/data_time.ascii")
# Timestep
dt = np.mean(time[1:] - time[:-1])

fig = plt.figure(figsize=(10, 7))
fig.set_tight_layout(True)


zoom_buffer = 0.028

rec_axwid = 0.35
rec_axhei = 0.3
z_hei = rec_axhei-0.1
z_wid = rec_axwid - 0.1
# Creating two axes
# add_axes([xmin,ymin,dx,dy])
ax_bl = fig.add_axes([0,0,rec_axwid,rec_axhei])
ax_br = fig.add_axes([1-rec_axwid,0,rec_axwid,rec_axhei])

ax_tl = fig.add_axes([0,1-rec_axhei,rec_axwid,rec_axhei])
ax_tr = fig.add_axes([1-rec_axwid,1-rec_axhei,rec_axwid,rec_axhei])

# Zoomed in:
ax_Zbl = fig.add_axes([0,rec_axhei-zoom_buffer,z_wid,z_hei])
ax_Zbr = fig.add_axes([1-z_wid,rec_axhei-zoom_buffer,z_wid,z_hei])

ax_Ztl = fig.add_axes([0,1-rec_axhei-z_hei + zoom_buffer,z_wid,z_hei])
ax_Ztr = fig.add_axes([1-z_wid,1-rec_axhei-z_hei + zoom_buffer,z_wid,z_hei])



# Clockwise list
ax_list = [ax_tr, ax_br, ax_bl, ax_tl]
Z_list  = [ax_Ztr, ax_Zbr, ax_Zbl, ax_Ztl]

for axi in [ax_bl, ax_br, ax_tl, ax_tr, ax_Ztr, ax_Ztl, ax_Zbl, ax_Zbr]:
    axi.spines[['right', 'top', 'left', 'bottom']].set_visible(False)
    # Remove axes ticks
    axi.tick_params(axis='both', which='both', bottom=False, left=False, top=False,
                labelbottom=False, labelleft=False)
    axi.patch.set_alpha(0)


# Add in the 3D model:
buffer = 0.1
arr_image = plt.imread('./Figures/model_clipped.png', format='png')
ax_model = fig.add_axes([z_wid-buffer, rec_axhei-buffer, 1-2*z_wid + 2*buffer, 1- 2*rec_axhei + 2*buffer], zorder=1)

ax_model.imshow(arr_image)
ax_model.spines[['right', 'top', 'left', 'bottom']].set_visible(False)
# Remove axes ticks
ax_model.tick_params(axis='both', which='both', bottom=False, left=False, top=False,
               labelbottom=False, labelleft=False)

ax_model.patch.set_alpha(0)




# Loop though each azimuth creating one figure for each azimuth
for azi in range(4):
    ax = ax_list[azi]
    axZ = Z_list[azi]
    for istn in range(1,nstns+1)[::-1]:
        # Produce obspy trace for filtering
        Ano = 2*azi + 1

        tr = Trace()
        tr.stats.delta = dt
        tr.data = np.loadtxt(f"{ddir}/{net}.A{Ano}_{istn}.ascii")[:,0]
        tr.filter('bandpass', freqmin=fmin, freqmax=fmax)

        # Convert trace to numpy arrays for mpl plotting
        t, d = obspy_gen_mpl(tr)

        # Some emirical scaling of data amplitude to separation on y axis
        if norm_traces:
            ax.plot(t, norm(d)+istn/1.75, color=hex[-istn], linewidth=0.8)
        else:
            ax.plot(t, d/0.08 +istn/1.75, color=hex[-istn], linewidth=0.8)

        if istn == 1 or istn==15:
            tarr = arr_time[istn-1]
            tmask = np.logical_and(t>=tarr-2.5, t<= tarr+20)
            axZ.plot(t[tmask]-tarr, 10*norm(d[tmask])+istn, color=hex[-istn], linewidth=1)
            axZ.set_xlim([-2.5, 17.5])

    # Some bells and whistles
    ax.set_ylim([-0.5, 10])

    ax.set_xlim([0, 50])
    #ax.set_title(f'Azimuth: {azi*45} degrees')
    #ax.spines[['right', 'top', 'left', 'bottom']].set_visible(False)
    # Remove axes ticks
    #ax.tick_params(axis='both', which='both', bottom=False, left=False, top=False,
    #                labelbottom=False, labelleft=False)

if save_figs_as_pdf:
    outstr = f'{odir}/bandpass_{fmin}_{fmax}_Hz'
    if norm_traces:
        outstr += '_norm'
    fig.savefig(f"{outstr}.pdf")
    print(f'Saved figure to {outstr}.pdf')

plt.show()
