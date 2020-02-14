import sys
from setup_figure import initialize_figure
from plotting_routines import plot_selfenergy, plot_spectralfunction

filename, cnt = initialize_figure()

plot = sys.argv[cnt]

fs = 14  # font size

if plot == 'selfenergy':
    plot_selfenergy(filename, cnt, fs, 13)
elif plot == 'spectralfunction':
    plot_spectralfunction(filename, fs, 0)

