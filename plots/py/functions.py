import os
import numpy as np
import yaml
from scipy.stats import iqr
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.cm as cm
import matplotlib.colors
import matplotlib.colorbar
import matplotlib.ticker
from mpl_toolkits.axes_grid1 import make_axes_locatable
from itertools import cycle
import time

# yaml
float_representer = lambda dumper,value: dumper.represent_scalar(u'tag:yaml.org,2002:float', "{:<.8e}".format(value))
npfloat_representer = lambda dumper,value: dumper.represent_float(float(value))
nparray_representer = lambda dumper,value: dumper.represent_list(value.tolist())
yaml.add_representer(float,float_representer)
yaml.add_representer(np.float_,npfloat_representer)
yaml.add_representer(np.ndarray,nparray_representer)

# global settings
plt.rcParams['axes.linewidth'] = 0.5

###########################################################################
## utils
###########################################################################
def print_time():
    return time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())

###########################################################################
## matrices
###########################################################################
def read_map(tfile,xy=False):
    # read matrices
    DATA=np.loadtxt(tfile)
    NN,D=DATA.shape

    if (D != 3):
        print "wrong format: D={:d}!".format(D)
        sys.exit()

    N = np.sqrt(float(NN))
    if (int(N) != N):
        print "Length is not a square root! sqrt(NN) = {:.4g}".format(N)
    N = int(N)

#    print "{:<20s}{:<d}".format("N",N)
#    print "{:<20s}{:<d}".format("D",D)

    MAT=np.reshape(DATA[:,2],(N,N))
    if not (xy):
        return MAT
    else:
        X=np.reshape(DATA[:,0],(N,N))
        Y=np.reshape(DATA[:,1],(N,N))
        return X,Y,MAT

def write_map(MAT,tfile,xy=True):
    if (len(MAT.shape) != 2):
        raise ValueError("wrong format: len(shape)={:d}!".format(len(MAT.shape)))

    N,N=MAT.shape

    fout = open(tfile,"w")
    for i in range(N):
        for j in range(N):
            val = MAT[i][j]
            if (xy):
                fout.write("{:<20d}{:<20d}{:<20.8e}\n".format(i,j,val))
            else:
                fout.write("{:<20.8e}\n".format(val))
    fout.close()
    return

def bin_map(x,y,mat,bin_size, binning='mean'):
    """
    return a binned version of the input matrix,
    with corresponding values for the x and y meshgrid.
    """
    bins = np.arange(np.min(x),np.max(x)+bin_size/2.0,bin_size)
    n = len(bins)-1
    digitized_x = np.digitize(x,bins)
    digitized_y = np.digitize(y,bins)
    mat_binned = np.zeros((n,n),dtype=np.float_)
    for i in range(n):
        for j in range(n):
            idx = ((digitized_x == i+1) & (digitized_y == j+1))
            if (binning == 'mean'):
                mat_binned[i,j] = np.mean(mat[idx])
            elif (binning == 'max'):
                mat_binned[i,j] = np.max(mat[idx])
            elif (binning == 'sum'):
                mat_binned[i,j] = np.sum(mat[idx])
    y,x = np.meshgrid(bins[:-1],bins[:-1])
    mat = mat_binned
    return x,y,mat

def compare_matrix(A,B,eucl_dist,rel_dist):
    if (len(A.shape) > 2):
        raise ValueError("A must be a matrix")
    if (A.shape != B.shape):
        raise ValueError("A and B must have same dimensions!")
    M = A.shape[0]
    N = A.shape[1]

    res = {}

    if (eucl_dist):
        res['eucl_dist'] = np.linalg.norm(A-B)/np.sqrt(N*M)
    if (rel_dist):
        res['rel_dist'] = (2.0*np.linalg.norm(A-B)/(np.linalg.norm(A)+np.linalg.norm(B)))
    return res

###########################################################################
## plots
###########################################################################

def plot_matrix(ax,X,Y,Z,cmapname,method=None,vmin=None,vmax=None, interpolation='none'):
    """
    Plot the input matrix into the given Axe instance.
    """
    if (vmin != None) and (vmax != None):
        print "vmin=%.6e  vmax=%.6f" %(vmin,vmax)
    if (cmapname == 'bone'):
        cmap=cm.bone_r
    if (cmapname == 'magma'):
        cmap=cm.magma_r
    if (cmapname == 'viridis'):
        cmap=cm.viridis_r
    elif (cmapname == 'jet'):
        cmap=cm.jet

    if (method == 'pcolormesh'):
        cs = ax.pcolormesh(X,Y,Z,cmap=cmap)
        xmin = np.min(X)
        xmax = np.max(X)
        print ax.set_xlim(xmin,xmax)
        print ax.set_ylim(xmin,xmax)
        ax.axis('equal')
        ax.invert_yaxis()
    else:
        cs = ax.imshow(Z,cmap=cmap, alpha=1.0, interpolation=interpolation, origin='upper', extent=[np.min(Y),np.max(Y),np.max(X),np.min(X)],vmin=vmin,vmax=vmax)
    return cs

def make_matrix_plot(filein, plot_matrix_arg, fig_width, fig_height, fontsize, nticks_max, bin_size,dpi,lognorm,lmin,amin,amax,arrows,circles,binning,cbarlabel=None):
    """
    Plot the matrix in the input file and write the output file.
    Argument:
        o plot_matrix_arg: dictionary of arguments to be passed to plot_matrix.
    """
    # import matrix
    x,y,mat = read_map(filein,xy=True)

    # make binning
    bin_size = int(bin_size)
    if (bin_size > 1):
        x,y,mat = bin_map(x,y,mat,bin_size=bin_size, binning=binning)


    # remove diagonal contacts
    idx = (np.abs(x-y) < lmin)
    mat[idx] = np.nan

    # remove extreme values
    if (amin != None):
        if (amin == 'median'):
            data = np.ravel(mat)
            idx = np.isfinite(data) & (data > 0.)
            data = data[idx]
            mat_med = np.median(data)
            mat_mad = 3*mad(data)
            amin = mat_med + 1*mat_mad
        elif (amin == 'nonzero'):
            data = np.ravel(mat)
            idx = np.isfinite(data) & (data > 0.)
            data = data[idx]
            amin = np.min(data)
#        amin = amin / bin_size**2
        idx = (mat < amin)
        mat[idx] = np.nan
    if (amax != None):
#        amax = amax / bin_size**2
        idx = (mat > amax)
        mat[idx] = np.nan

    # create figure
    fig=plt.figure(num='none',facecolor='white', figsize=(fig_width,fig_height))
    ax=fig.gca()

    # plot figure
    cs=plot_matrix(ax,x,y,mat,**plot_matrix_arg)

    # draw arrows on largest values
    if (arrows != None):
        narrows = arrows['n']
        angle_annot = arrows['theta'] * np.pi / 180.0
        r_annot = arrows['l']
        xflat = np.ravel(x)
        yflat = np.ravel(y)
        mflat = np.ravel(mat)
        idx = np.isfinite(mflat)
        xflat = xflat[idx]
        yflat = yflat[idx]
        mflat = mflat[idx]
        if ('amin' in arrows) and (arrows['amin'] != None):
            mmin = arrows['amin']
            mmin = mmin / bin_size**2
            idx = (mflat >= mmin)
            xflat = xflat[idx]
            yflat = yflat[idx]
            mflat = mflat[idx]
        if ('amax' in arrows) and (arrows['amax'] != None):
            mmax = arrows['amax']
            mmax = mmax / bin_size**2
            idx = (mflat <= mmax)
            xflat = xflat[idx]
            yflat = yflat[idx]
            mflat = mflat[idx]
        if ('lmin' in arrows) and (arrows['lmin'] != None):
            idx = (np.abs(xflat-yflat) >= arrows['lmin'])
            xflat = xflat[idx]
            yflat = yflat[idx]
            mflat = mflat[idx]
        # select only one side
        idx = (xflat-yflat > 0)
        xflat = xflat[idx]
        yflat = yflat[idx]
        mflat = mflat[idx]
        idx = np.argsort(mflat)[::-1]
        xflat = xflat[idx]
        yflat = yflat[idx]
        mflat = mflat[idx]
        narrows=min(len(mflat),narrows)
        for na in range(narrows):
            x_annot = xflat[na]
            y_annot = yflat[na]
            xy_annot = (x_annot,y_annot)
            xy_text = (r_annot*np.cos(angle_annot),r_annot*np.sin(angle_annot))

            ax.annotate('', xy=xy_annot, xycoords='data',
                            xytext = xy_text, textcoords='offset points',
                        arrowprops=dict(arrowstyle='-|>',facecolor='black', shrinkA=0.00, shrinkB=2.00,lw=0.5)
                        )
    if (circles != None):
        ncircles = circles['n']
        r_annot = circles['r']
        xflat = np.ravel(x)
        yflat = np.ravel(y)
        mflat = np.ravel(mat)
        idx = np.isfinite(mflat)
        xflat = xflat[idx]
        yflat = yflat[idx]
        mflat = mflat[idx]
        if ('amin' in circles) and (circles['amin'] != None):
            mmin = circles['amin']
            mmin = mmin / bin_size**2
            idx = (mflat >= mmin)
            xflat = xflat[idx]
            yflat = yflat[idx]
            mflat = mflat[idx]
        if ('amax' in circles) and (circles['amax'] != None):
            mmax = circles['amax']
            mmax = mmax / bin_size**2
            idx = (mflat <= mmax)
            xflat = xflat[idx]
            yflat = yflat[idx]
            mflat = mflat[idx]
        if ('lmin' in circles) and (circles['lmin'] != None):
            idx = (np.abs(xflat-yflat) >= arrows['lmin'])
            idx = (np.abs(xflat-yflat) >= circles['lmin'])
            xflat = xflat[idx]
            yflat = yflat[idx]
            mflat = mflat[idx]
        # select only one side
        idx = (xflat-yflat > 0)
        xflat = xflat[idx]
        yflat = yflat[idx]
        mflat = mflat[idx]
        idx = np.argsort(mflat)[::-1]
        xflat = xflat[idx]
        yflat = yflat[idx]
        mflat = mflat[idx]
        ncircles=min(len(mflat),ncircles)
        for nc in range(ncircles):
            x_annot = xflat[nc]
            y_annot = yflat[nc]
            xy_annot = (x_annot,y_annot)
            circle = plt.Circle(xy_annot, r_annot, color='k',fill=False,lw=0.5)
            ax.add_artist(circle)

    # customize axis
    tick_locator = matplotlib.ticker.MaxNLocator(nbins=nticks_max,integer=True)
    ax.xaxis.set_major_locator(tick_locator)
    ax.yaxis.set_major_locator(tick_locator)
    ax.tick_params(length=2)

    # customize color bar
    if (lognorm):
        cs.norm = matplotlib.colors.LogNorm()
        tick_locator = matplotlib.ticker.LogLocator(subs=np.arange(10,dtype=np.float_)/10.,numticks=8)
        tick_formatter = matplotlib.ticker.LogFormatterMathtext()
    else:
        cs.norm = matplotlib.colors.Normalize()
        tick_locator = matplotlib.ticker.AutoLocator()
        tick_formatter = matplotlib.ticker.ScalarFormatter()

    divider = make_axes_locatable(ax)
    cbar_ax = divider.append_axes("right", "5%", pad="10%")
    cbar=fig.colorbar(cs,cax=cbar_ax,ticks=tick_locator,format=tick_formatter)
    cbar.ax.tick_params(length=2)
    if (cbarlabel != None):
        cbar.ax.set_ylabel(cbarlabel,fontsize=fontsize, rotation=270, labelpad=20)

    # legend
    rect=[0.0,0.0,1.0,1.0]
    fig.tight_layout(rect=rect,pad=0.1)
    destdir=os.path.dirname(os.path.relpath(filein))
    basename = os.path.splitext(os.path.basename(filein))[0]
    fileout=os.path.join(destdir,"{}".format(basename))
    if (bin_size > 1):
        fileout = "{}_bin{:d}".format(fileout,bin_size)
#    if (lmin > 0):
#        fileout = "{}_lmin{:d}".format(fileout,lmin)
    if (lognorm):
        fileout="{}_lognorm".format(fileout)
#    fileout = "{}.png".format(fileout)
    fileout = "{}.pdf".format(fileout)
    fig.savefig(fileout,dpi=dpi,bbox_inches='tight',pad_inches=0)
    plt.close('all')

    return fileout

def make_matrix_plot_two(file_exp, file_pred, plot_matrix_arg, fig_width, fig_height, fontsize, nticks_max, bin_size,dpi,lognorm,lmin,amin,amax,arrows,circles,binning):
    """
    Plot the matrix in the input file and write the output file.
    Argument:
        o plot_matrix_arg: dictionary of arguments to be passed to plot_matrix.
    """
    # import matrix
    x,y,mat_exp = read_map(file_exp,xy=True)
    x,y,mat_pred = read_map(file_pred,xy=True)

    # make binning
    bin_size = int(bin_size)
    if (bin_size > 1):
        x,y,mat_exp = bin_map(x,y,mat_exp,bin_size=bin_size,binning=binning)
        x,y,mat_pred = bin_map(x,y,mat_pred,bin_size=bin_size,binning=binning)

    # combine matrices
    tri_lo = (x >= y)
    tri_up = (x < y)
    mat = np.empty(mat_exp.shape, dtype=np.float_)
    mat[tri_lo] = mat_exp[tri_lo]
    mat[tri_up] = mat_pred[tri_up]

    # remove diagonal contacts
    idx = (np.abs(x-y) < lmin)
    mat[idx] = np.nan

    # remove extreme values
    if (amin != None):
        if (amin == 'median'):
            data = np.ravel(mat)
            idx = np.isfinite(data) & (data > 0.)
            data = data[idx]
            mat_med = np.median(data)
            mat_mad = 3*mad(data)
            amin = mat_med + 1*mat_mad
        elif (amin == 'nonzero'):
            data = np.ravel(mat)
            idx = np.isfinite(data) & (data > 0.)
            data = data[idx]
            amin = np.min(data)
#        amin = amin / bin_size**2
        idx = (mat < amin)
        mat[idx] = np.nan
    if (amax != None):
#        amax = amax / bin_size**2
        idx = (mat > amax)
        mat[idx] = np.nan

    # create figure
    fig=plt.figure(num='none',facecolor='white', figsize=(fig_width,fig_height))
    ax=fig.gca()

    # plot figure
    cs=plot_matrix(ax,x,y,mat,**plot_matrix_arg)

    # draw arrows on largest values
    if (arrows != None):
        narrows = arrows['n']
        angle_annot = arrows['theta'] * np.pi / 180.0
        r_annot = arrows['l']
        xflat = np.ravel(x)
        yflat = np.ravel(y)
        mflat = np.ravel(mat)
        idx = np.isfinite(mflat)
        xflat = xflat[idx]
        yflat = yflat[idx]
        mflat = mflat[idx]
        if ('amin' in arrows) and (arrows['amin'] != None):
            mmin = arrows['amin']
            mmin = mmin / bin_size**2
            idx = (mflat >= mmin)
            xflat = xflat[idx]
            yflat = yflat[idx]
            mflat = mflat[idx]
        if ('amax' in arrows) and (arrows['amax'] != None):
            mmax = arrows['amax']
            mmax = mmax / bin_size**2
            idx = (mflat <= mmax)
            xflat = xflat[idx]
            yflat = yflat[idx]
            mflat = mflat[idx]
        if ('lmin' in arrows) and (arrows['lmin'] != None):
            idx = (np.abs(xflat-yflat) >= arrows['lmin'])
            xflat = xflat[idx]
            yflat = yflat[idx]
            mflat = mflat[idx]
        # select only one side
        idx = (xflat-yflat > 0)
        xflat = xflat[idx]
        yflat = yflat[idx]
        mflat = mflat[idx]
        idx = np.argsort(mflat)[::-1]
        xflat = xflat[idx]
        yflat = yflat[idx]
        mflat = mflat[idx]
        narrows=min(len(mflat),narrows)
        for na in range(narrows):
            x_annot = xflat[na]
            y_annot = yflat[na]
            xy_annot = (x_annot,y_annot)
            xy_text = (r_annot*np.cos(angle_annot),r_annot*np.sin(angle_annot))

            ax.annotate('', xy=xy_annot, xycoords='data',
                            xytext = xy_text, textcoords='offset points',
                        arrowprops=dict(arrowstyle='-|>',facecolor='black', shrinkA=0.00, shrinkB=2.00,lw=0.5)
                        )
    if (circles != None):
        ncircles = circles['n']
        r_annot = circles['r']
        xflat = np.ravel(x)
        yflat = np.ravel(y)
        mflat = np.ravel(mat)
        idx = np.isfinite(mflat)
        xflat = xflat[idx]
        yflat = yflat[idx]
        mflat = mflat[idx]
        if ('amin' in circles) and (circles['amin'] != None):
            mmin = circles['amin']
            mmin = mmin / bin_size**2
            idx = (mflat >= mmin)
            xflat = xflat[idx]
            yflat = yflat[idx]
            mflat = mflat[idx]
        if ('amax' in circles) and (circles['amax'] != None):
            mmax = circles['amax']
            mmax = mmax / bin_size**2
            idx = (mflat <= mmax)
            xflat = xflat[idx]
            yflat = yflat[idx]
            mflat = mflat[idx]
        if ('lmin' in circles) and (circles['lmin'] != None):
            idx = (np.abs(xflat-yflat) >= arrows['lmin'])
            idx = (np.abs(xflat-yflat) >= circles['lmin'])
            xflat = xflat[idx]
            yflat = yflat[idx]
            mflat = mflat[idx]
        # select only one side
        idx = (xflat-yflat > 0)
        xflat = xflat[idx]
        yflat = yflat[idx]
        mflat = mflat[idx]
        idx = np.argsort(mflat)[::-1]
        xflat = xflat[idx]
        yflat = yflat[idx]
        mflat = mflat[idx]
        ncircles=min(len(mflat),ncircles)
        for nc in range(ncircles):
            x_annot = xflat[nc]
            y_annot = yflat[nc]
            xy_annot = (x_annot,y_annot)
            circle = plt.Circle(xy_annot, r_annot, color='k',fill=False,lw=0.5)
            ax.add_artist(circle)

    # customize axis
    tick_locator = matplotlib.ticker.MaxNLocator(nbins=nticks_max,integer=True)
    ax.xaxis.set_major_locator(tick_locator)
    ax.yaxis.set_major_locator(tick_locator)
    ax.tick_params(length=2)

    # customize color bar
    if (lognorm):
        cs.norm = matplotlib.colors.LogNorm()
        tick_locator = matplotlib.ticker.LogLocator(subs=np.arange(10,dtype=np.float_)/10.,numticks=8)
        tick_formatter = matplotlib.ticker.LogFormatterMathtext()
    else:
        cs.norm = matplotlib.colors.Normalize()
        tick_locator = matplotlib.ticker.AutoLocator()
        tick_formatter = matplotlib.ticker.ScalarFormatter()

    divider = make_axes_locatable(ax)
    cbar_ax = divider.append_axes("right", "5%", pad="10%")
    cbar=fig.colorbar(cs,cax=cbar_ax,ticks=tick_locator,format=tick_formatter)

    # legend
    rect=[0.0,0.0,1.0,1.0]
    fig.tight_layout(rect=rect,pad=0.1)
    destdir=os.path.dirname(os.path.relpath(file_pred))
    basename_exp = os.path.splitext(os.path.basename(file_exp))[0]
    basename_pred = os.path.splitext(os.path.basename(file_pred))[0]
    fileout=os.path.join(destdir,"{}_{}".format(basename_exp,basename_pred))
    if (bin_size > 1):
        fileout = "{}_bin{:d}".format(fileout,bin_size)
    if (lmin > 0):
        fileout = "{}_lmin{:d}".format(fileout,lmin)
    if (lognorm):
        fileout="{}_lognorm".format(fileout)
#    fileout = "{}.png".format(fileout)
    fileout = "{}.pdf".format(fileout)
    fig.savefig(fileout,dpi=dpi,bbox_inches='tight',pad_inches=0)
    plt.close('all')

    return fileout

def make_matrix_plot_multi(fileins, plot_matrix_arg, ax_width, ax_height, fontsize, nticks_max, bin_size, dpi, lognorm, nrows, ncols, amin, amax, titles=[], ref=None, min_sub='none', max_sub='none', cmat_ref_file=None):
    """
    Plot the matrix in the input file and write the output file.
    Argument:
        o plot_matrix_arg: dictionary of arguments to be passed to plot_matrix.

    For colorbars and norms, see doc at:
        o http://matplotlib.org/examples/api/colorbar_only.html
        o http://matplotlib.org/api/cm_api.html
    """
    data=[]
    basenames = []
    vmin=99.9e99
    vmax=-99.9e99
    mat_ref = None
    if cmat_ref_file != None:
        x,y,mat_ref = read_map(cmat_ref_file,xy=True)

    for filein in fileins:
        # import matrices
        x,y,mat_raw = read_map(filein,xy=True)

        # combine matrices
        if cmat_ref_file != None:
            tri_lo = (x >= y)
            tri_up = (x < y)
            mat = np.empty(mat_ref.shape, dtype=np.float_)
            mat[tri_lo] = mat_ref[tri_lo]
            mat[tri_up] = mat_raw[tri_up]

        # make binning
        bin_size = int(bin_size)
        if (bin_size > 1):
            x,y,mat = bin_map(x,y,mat,bin_size=bin_size,binning='mean')
        else:
            mat = mat_raw

        # remove extreme values
        if (amin != None):
#            amin = amin / bin_size**2
            idx = np.logical_not(np.isnan(mat)) & (mat < amin)
            if (min_sub == 'none'):
                mat[idx]=np.nan
            else:
                mat[idx] = amin
        if (amax != None):
#            amax = amax / bin_size**2
            idx = np.logical_not(np.isnan(mat)) & (mat > amax)
            if (max_sub == 'none'):
                mat[idx]=np.nan
            else:
                mat[idx] = vmax

        vmin = min(np.min(mat[np.isfinite(mat)]),vmin)
        vmax = max(np.max(mat[np.isfinite(mat)]),vmax)

        # add matrix
        data.append((x,y,mat))

        # basenames
        basenames.append(os.path.splitext(os.path.basename(filein))[0])

    # create figure
    fig=plt.figure(num='none',facecolor='white', figsize=(ncols*ax_width,nrows*ax_height))
    gs = gridspec.GridSpec(nrows,ncols)
    gs.update(left=0.00, right=0.95, wspace=0.0)
    gs_cbar = gridspec.GridSpec(1,1)
    gs_cbar.update(left=0.98, right=1.00, wspace=0)
    css = []

    # choose reference plot
    selection=range(len(data))
    if (ref != None):
        selection=range(ref)+range(ref+1,len(data))
        x,y,matref=data[ref]
        idx=(x>y)
        for i in selection:
            x,y,mat=data[i]
            mat[idx]=matref[idx]
            data[i]=x,y,mat

    # loop on plots
    for i,ind in enumerate(selection):
        # control
        if not (i < nrows*ncols ):
            print "There are not enough plots in the figures!"
            break

        # get data
        x,y,mat = data[ind]

        # plot figure
        c = i%ncols
        r = (i-c)/ncols
        ax = fig.add_subplot(gs[r,c])
        cs=plot_matrix(ax,x,y,mat,vmin=vmin,vmax=vmax,**plot_matrix_arg)
        css.append(cs)

        if (lognorm):
            cs.norm = matplotlib.colors.LogNorm(vmin=vmin,vmax=vmax)
        else:
            cs.norm = matplotlib.colors.Normalize(vmin=vmin,vmax=vmax)

        # customize axis
        tick_locator_x = matplotlib.ticker.MaxNLocator(nbins=nticks_max,integer=True)
        ax.xaxis.set_major_locator(tick_locator_x)
        ax.yaxis.set_major_locator(tick_locator_x)
        ax.tick_params(length=2)

        # title
        if (titles != []):
            text = titles[ind]
            xy = (0.5,1.0)
            ax.annotate(text,xy=xy,xycoords='axes fraction',ha='center',va='bottom', fontsize=1.2*fontsize)

    # customize color bar
    cbar_ax = fig.add_subplot(gs_cbar[0,0])
    if (lognorm):
        tick_locator = matplotlib.ticker.LogLocator(subs=np.arange(10,dtype=np.float_)/10.,numticks=8)
        tick_formatter = matplotlib.ticker.LogFormatterMathtext()
    else:
        tick_locator = matplotlib.ticker.AutoLocator()
        tick_formatter = matplotlib.ticker.ScalarFormatter()

    cbar = matplotlib.colorbar.Colorbar(ax=cbar_ax, mappable=cs, ticks=tick_locator, format=tick_formatter)

    # legend
    destdir=os.path.dirname(os.path.relpath(fileins[0]))
    suf = "_".join(basenames)
    fileout=os.path.join(destdir,"{}".format(suf))
    if (bin_size > 1):
        fileout = "{}_bin{:d}".format(fileout,bin_size)
    if (lognorm):
        fileout="%s_lognorm" %(fileout)
#    fileout = "{}.png".format(fileout)
    fileout = "{}.pdf".format(fileout)
    fig.savefig(fileout,dpi=dpi,bbox_inches='tight',pad_inches=0)
    plt.close('all')

    return fileout

def make_matrix_compare_multi(fileins, compare_matrix_arg):
    """
    Plot the matrix in the input file and write the output file.
    Argument:
        o plot_matrix_arg: dictionary of arguments to be passed to plot_matrix.

    For colorbars and norms, see doc at:
        o http://matplotlib.org/examples/api/colorbar_only.html
        o http://matplotlib.org/api/cm_api.html
    """

    if (fileins == []):
        return None

    res = {}
    basenames=[]
    # loop and import matrices
    nf = len(fileins)
    for i in range(nf):
        filein1 = fileins[i]
        # basenames
        basename1 = os.path.splitext(os.path.basename(filein1))[0]
        basenames.append(basename1)

        if (i == nf -1):
            break

        res[basename1]={}

        # read data
        mat1 = read_map(filein1)

        for j in range(i+1,nf):
            filein2 = fileins[j]
            basename2 = os.path.splitext(os.path.basename(filein2))[0]
            mat2 = read_map(filein2)
            res[basename1][basename2] = compare_matrix(mat1,mat2,**compare_matrix_arg)

    # write result and return
    destdir=os.path.dirname(os.path.relpath(fileins[0]))
    suf = "_".join(basenames)
    fileout=os.path.join(destdir,"{}".format(suf))
    fileout = "{}.dat".format(fileout)
    fout = open(fileout,"w")
    yaml.dump(res,stream=fout,default_flow_style=False, tags=None)
    fout.close()
    return fileout
