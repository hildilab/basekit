from __future__ import division


import math

from utils import listify, flatten, get_index

import numpy as np

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import cm


def rama_plot( phi_psi, suptitle="Ramachandran Plot", titles=None,
               xlabel="Psi", ylabel="Phi", image_file="rama.png" ):
    """
        phi_psi: 
            [ (phi1, psi1), (phi2, psi2), ... ] 
            or (phi1, psi1)
        titles: 
            [ title1, title2, ... ] 
            or title1
    """
    phi_psi = listify( phi_psi )
    titles = listify( titles )
    n = len( phi_psi )
    if not titles and n>1:
        titles = map( str, range(n) )

    nrows = int( round( math.sqrt(n) ) )
    ncols = int( math.ceil( math.sqrt(n) ) )
    figsize = ( ncols*5, nrows*5 )

    fig, axes = plt.subplots( nrows=nrows, ncols=ncols, figsize=figsize)
    fig.suptitle( suptitle )

    for i in range(n):
        if n==1:
            ax = axes
        else:
            ax = axes[ int( math.floor( i/ncols ) ), i%ncols ]
        phi = np.ma.masked_invalid( phi_psi[i][0], copy=False )
        psi = np.ma.masked_invalid( phi_psi[i][1], copy=False )
        _rama_sub_plot( 
            phi, psi, ax, xlabel=xlabel, ylabel=ylabel, 
            title=get_index( titles, i, "" )
        )
    
    fig.tight_layout()
    fig.savefig( image_file, dpi=300 )


def _rama_sub_plot( phi, psi, ax, xlabel="Psi", ylabel="Phi", title=None ):
    xy_range = [ [ -180, 180 ], [ -180, 180 ] ]
    hist2d, phi_edges, psi_edges = np.histogram2d( 
        phi, psi, bins=180, range=xy_range
    )
    ax.imshow( 
        np.flipud( hist2d ), extent=flatten( xy_range ), 
        interpolation='nearest', cmap=cm.Reds
    )
    ax.set_title( title )
    ax.set_xlabel( xlabel )
    ax.set_ylabel( ylabel )
