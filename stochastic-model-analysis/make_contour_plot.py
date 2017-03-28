import matplotlib
import numpy as np
import matplotlib.cm as cm
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import scipy.interpolate

font = {'family': 'Arial',
        'weight': 'bold',
        'size': 10,
        }
####LOAD DATA ####

obs_data_full = np.genfromtxt('model-data/contour_plot_data.csv', delimiter=',', skip_header=1)

####Interpolation for contours and heatmap
R0 = obs_data_full[:,0]
I0 = np.log10(obs_data_full[:,1])
prob = obs_data_full[:,4]

xi_full, yi_full2 = np.linspace(I0.min(), 2, 300), np.linspace(R0.min(), 1.5, 300)
xi_full, yi_full2 = np.meshgrid(xi_full, yi_full2)
zi_full2 = scipy.interpolate.griddata((I0, R0), prob, (xi_full, yi_full2), method='linear')


#### Make plot
fig3 = plt.subplot(1,1,1)
plt.contourf(xi_full, yi_full2, zi_full2, 9, cmap=plt.cm.jet, alpha = 0.75)
plt.colorbar() 
CS = plt.contour(xi_full, yi_full2, zi_full2, 9, colors='black', linewidth=.5)


# for labels but plotted on the log scale
class nf(float):
    def __repr__(self):
        str = '%.1f' % (self.__float__(),)
        if str[-1] == '0':
            return '%.0f' % self.__float__()
        else:
            return '%.1f' % self.__float__()

# Recast levels to new class
CS.levels = [nf(val) for val in CS.levels]
# Label levels with specially formatted floats
if plt.rcParams["text.usetex"]:
    fmt = r'%r \%%'
else:
    fmt = '%r'

plt.clabel(CS, CS.levels, inline=True, inline_spacing = 12, fontsize=12, fmt = fmt, color = 'black')

a=fig3.get_xticks().tolist()
for i in [0,2,4]:
	a[i] = 10**(i/2)
for i in [1,3]: 
	a[i] = ''
fig3.set_xticklabels(a)

### Added dotted line at 0.76
CS = plt.contour(xi_full, yi_full2, zi_full2, levels = [0.76], colors='black', linewidth=.5, linestyles = "dotted")


plt.xlabel(r'$I_0$', fontsize=14)
plt.ylabel(r'$R_0$', fontsize=14)
fig3 = matplotlib.pyplot.gcf()
plt.tight_layout()
plt.savefig('transient_infection_prob_heatmap.pdf')

plt.close()
