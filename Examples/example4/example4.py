from MOOSEplot import geoplot as gp
from MOOSEplot import plotCommon as pc
from MOOSEfunctions import ColormapLoad
from matplotlib import pyplot as plt

cmaps=ColormapLoad.load_all_cmaps(verbose=True)
Gplot=gp.StandardPlot()
Gplot.parse_namelist(namelist_file='radar.txt')
Gplot.namelist_dictionary['colorbar']=False

Gplot.define_projection()
fig=plt.figure(figsize=(8,12.5))
ax=plt.subplot(1,1,1,projection=Gplot.projection)

Gplot.makePlot(ax)

pos=ax.get_position()
cbar_ax=fig.add_axes([pos.x0+0.03,pos.y1-0.05,0.6,0.017])
cbar=plt.colorbar(Gplot.im, cax=cbar_ax,ticks=[-10,0,10,20,30,40,50,60,70,80],orientation='horizontal')
cbar.set_label("dBZ")
plt.show()
