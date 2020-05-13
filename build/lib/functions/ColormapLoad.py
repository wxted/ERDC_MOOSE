# coding=utf-8
"""This file contains a few helper functions to generate new color
maps based on input text (.rgb or .txt.) files with RGB color information.
Admittedly, there are not enough comments or documentation here, and a lot of the code was more or less
copied over from Stack Overflow."""

from matplotlib import pyplot as plt
import numpy as np
from matplotlib import colors as mplcol
from matplotlib import cm as mplcm
import glob as glob
import os

def cmap_discretize(cmap, N):
    """Return a discrete colormap from the continuous colormap cmap.

        cmap: colormap instance, eg. cm.jet.
        N: number of colors.

    Example
        x = resize(arange(100), (5,100))
        djet = cmap_discretize(cm.jet, 5)
        imshow(x, cmap=djet)
    """
    if type(cmap) == str:
        cmap = plt.get_cmap(cmap)
    colors_i = np.concatenate((np.linspace(0, 1., N), (0.,0.,0.,0.)))
    colors_rgba = cmap(colors_i)
    indices = np.linspace(0, 1., N+1)
    cdict = {}

    print(len(indices))
    for ki,key in enumerate(('red','green','blue')):
        cdict[key] = [ (indices[i], colors_rgba[i-1,ki], colors_rgba[i,ki])
                       for i in range(N+1) ]
 # Return colormap object.
    return mplcol.LinearSegmentedColormap(cmap.name + "_%d"%N, cdict, 1024)



def mpy_to_mpl_cmap(ctable,rev=False):
        """This function grabs a colortable with the metpy format and converts it to a table that can be easily parsed as an mpl map"""

        LinL=[line.replace('(', '').replace(')', '').rstrip('\n') for line in open(ctable)]
        LinL=[np.fromstring(l, dtype=float, sep=",") for l in LinL]
        LinL=np.vstack(LinL)

        if rev == True:
            LinL=LinL[::-1]

        b3=LinL[:,2] # value of blue at sample n
        b2=LinL[:,2] # value of blue at sample n
        b1=np.linspace(0,1,len(b2)) # position of sample n - ranges from 0 to 1

        # setting up columns for list
        g3=LinL[:,1]
        g2=LinL[:,1]
        g1=np.linspace(0,1,len(g2))

        r3=LinL[:,0]
        r2=LinL[:,0]
        r1=np.linspace(0,1,len(r2))

        # creating list
        R=zip(r1,r2,r3)
        G=zip(g1,g2,g3)
        B=zip(b1,b2,b3)

        # transposing list
        RGB=zip(R,G,B)
        rgb=zip(*RGB)
        # print rgb
        # creating dictionary
        k=['red', 'green', 'blue']
        LinearL=dict(zip(k,rgb)) # makes a dictionary from 2 lists
        my_cmap = mplcol.LinearSegmentedColormap('my_colormap',LinearL)
        return my_cmap



def ncl_to_mpl_cmap(ctable,rev=False):
#        path='colortables/'
        #__location__=os.path.dirname(os.path.abspath(__file__))
       # path = os.path.join(__location__, path)
        LinL=np.loadtxt(ctable,comments='#',skiprows=2)/255.

        if rev == True:
            LinL=LinL[::-1]

        b3=LinL[:,2] # value of blue at sample n
        b2=LinL[:,2] # value of blue at sample n
        b1=np.linspace(0,1,len(b2)) # position of sample n - ranges from 0 to 1

        # setting up columns for list
        g3=LinL[:,1]
        g2=LinL[:,1]
        g1=np.linspace(0,1,len(g2))

        r3=LinL[:,0]
        r2=LinL[:,0]
        r1=np.linspace(0,1,len(r2))

        # creating list
        R=zip(r1,r2,r3)
        G=zip(g1,g2,g3)
        B=zip(b1,b2,b3)

        # transposing list
        RGB=zip(R,G,B)
        rgb=zip(*RGB)
        # print rgb
        # creating dictionary
        k=['red', 'green', 'blue']
        LinearL=dict(zip(k,rgb)) # makes a dictionary from 2 lists
        my_cmap = mplcol.LinearSegmentedColormap('my_colormap',LinearL)
        return my_cmap


def load_all_cmaps(verbose=True,pp='Auto'):
    path=os.getcwd()
    #__location__=os.path.dirname(os.path.abspath(__file__))
#    path = os.path.join(__location__, path)
    print("LOADING COLORMAPS!")
    if pp.lower() == 'auto':
        cmaps=glob.glob(path+'colortables/*.rgb')
    else:
        cmaps=glob.glob(pp+'/*.rgb')
    for m in cmaps:
        name=m.split('/')[-1].split('.rgb')[0]
        if verbose == True:
            print ("Colormap: "+name+" loaded.")
        my_cmap=ncl_to_mpl_cmap(m)
        mplcm.register_cmap(name=name, cmap=my_cmap)
        my_cmap=ncl_to_mpl_cmap(m,rev=True)
        mplcm.register_cmap(name=name+'_r', cmap=my_cmap)

    ## NOW CONVERT OTHER TABLES! ##
    if pp.lower() == 'auto':
        cmaps=glob.glob(path+'colortables/*.tbl')
    else:
        cmaps=glob.glob(pp+'/*.tbl')
    for m in cmaps:
        name=m.split('/')[-1].split('.tbl')[0]
        if verbose == True:
            print ("Colormap: "+name+" loaded.")
        my_cmap=mpy_to_mpl_cmap(m)
        mplcm.register_cmap(name=name, cmap=my_cmap)
        my_cmap=mpy_to_mpl_cmap(m,rev=True)
        mplcm.register_cmap(name=name+'_r', cmap=my_cmap)
