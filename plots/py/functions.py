#################### color list ####################
# indigo, gold, hotpink, firebrick, indianred, yellow, mistyrose, darkolivegreen, olive, darkseagreen, pink, tomato, lightcoral, orangered, navajowhite, lime, palegreen, darkslategrey, greenyellow, burlywood, seashell, mediumspringgreen, fuchsia, papayawhip, blanchedalmond, chartreuse, dimgray, black, peachpuff, springgreen, aquamarine, white, orange, lightsalmon, darkslategray, brown, ivory, dodgerblue, peru, lawngreen, chocolate, crimson, forestgreen, darkgrey, lightseagreen, cyan, mintcream, silver, antiquewhite, mediumorchid, skyblue, gray, darkturquoise, goldenrod, darkgreen, floralwhite, darkviolet, darkgray, moccasin, saddlebrown, grey, darkslateblue, lightskyblue, lightpink, mediumvioletred, slategrey, red, deeppink, limegreen, darkmagenta, palegoldenrod, plum, turquoise, lightgrey, lightgoldenrodyellow, darkgoldenrod, lavender, maroon, yellowgreen, sandybrown, thistle, violet, navy, magenta, dimgrey, tan, rosybrown, olivedrab, blue, lightblue, ghostwhite, honeydew, cornflowerblue, slateblue, linen, darkblue, powderblue, seagreen, darkkhaki, snow, sienna, mediumblue, royalblue, lightcyan, green, mediumpurple, midnightblue, cornsilk, paleturquoise, bisque, slategray, darkcyan, khaki, wheat, teal, darkorchid, salmon, deepskyblue, rebeccapurple, darkred, steelblue, palevioletred, lightslategray, aliceblue, lightslategrey, lightgreen, orchid, gainsboro, mediumseagreen, lightgray, mediumturquoise, lemonchiffon, cadetblue, lightyellow, lavenderblush, coral, purple, aqua, whitesmoke, mediumslateblue, darkorange, mediumaquamarine, darksalmon, beige, blueviolet, azure, lightsteelblue, oldlace
#################### imports ####################
import os,sys
import numpy as np
import yaml
from scipy.stats import iqr,pearsonr
import scipy.stats as sstat
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

# matplotlib controls
plt.rcParams['svg.fonttype'] = 'none'  # to embed fonts in output ('path' is to convert as text as paths)
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams['axes.linewidth']=0.5
colornames = matplotlib.colors.cnames

# dpi
mydpi=150

#################### functions -- utils ####################
def print_time():
    return time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())

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

#################### functions -- plots ####################
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
    cs=plot_matrix(ax,x,y,mat,interpolation='nearest',**plot_matrix_arg) # note that interpolation='none' causes a problem with the svg image type

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
        cbar.ax.set_ylabel(cbarlabel,fontsize='large', rotation=270, labelpad=20)

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
    filename = fileout
    for ext in '.png', '.pdf', '.svg':
        fileout = filename + ext
        fig.savefig(fileout,bbox_inches='tight',pad_inches=0,dpi=mydpi)
        print "{:<20s}{:<s}".format("fileout",fileout)
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
#        mat[idx] = amin
    if (amax != None):
#        amax = amax / bin_size**2
        idx = (mat > amax)
        mat[idx] = np.nan

    # create figure
    fig=plt.figure(num='none',facecolor='white', figsize=(fig_width,fig_height))
    ax=fig.gca()

    # plot figure
    cs=plot_matrix(ax,x,y,mat,interpolation='nearest',**plot_matrix_arg) # note that interpolation='none' causes a problem with the svg image type

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
    filename=fileout
    for ext in '.png', '.pdf', '.svg':
        fileout = filename + ext
        fig.savefig(fileout,bbox_inches='tight',pad_inches=0,dpi=mydpi)
        print "{:<20s}{:<s}".format("fileout",fileout)
    plt.close('all')

    return fileout

def make_matrix_plot_multi(fileins, plot_matrix_arg, ax_width, ax_height, fontsize, nticks_max, bin_size, dpi, lognorm, nrows, ncols, amin, amax, titles=[], ref=None, min_sub='none', max_sub='none', cmat_ref_file=None,filename=None):
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
    gs.update(left=0.00, right=0.98, wspace=0.0)
    gs_cbar = gridspec.GridSpec(1,1)
    gs_cbar.update(left=0.985, right=1.00, wspace=0)
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
        print "i = {:d}    ind = {:d}".format(i,ind)
        # control
        if not (i < nrows*ncols ):
            print "There are not enough axes in the figures!"
            break

        # get data
        x,y,mat = data[ind]

        # plot figure
        c = i%ncols
        r = (i-c)/ncols
        ax = fig.add_subplot(gs[r,c])
        cs=plot_matrix(ax,x,y,mat,vmin=vmin,vmax=vmax,interpolation='nearest',**plot_matrix_arg)
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
            ax.annotate(text,xy=xy,xycoords='axes fraction',ha='center',va='bottom', fontsize='medium')

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
    if filename is None:
        destdir=os.path.dirname(os.path.relpath(fileins[0]))
        suf = "_".join(basenames)
        filename=os.path.join(destdir,"{}".format(suf))
    if (bin_size > 1):
        filename = "{}_bin{:d}".format(filename,bin_size)
    if (lognorm):
        filename="%s_lognorm" %(filename)
    for ext in '.png', '.pdf', '.svg':
        fileout = filename + ext
        fig.savefig(fileout,bbox_inches='tight',pad_inches=0,dpi=mydpi)
        print "{:<20s}{:<s}".format("fileout",fileout)
    #fig.savefig(fileout,dpi=dpi,transparent=True, bbox_inches='tight',pad_inches=0)
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

def make_distance_plot(filein, col_thres, col_dist1, col_dist2, xlabel, ylabel1, ylabel2, annot=True, ticks_dx=0.5, ticks_dy=0.1, lw=0.5, ms=2, linestyle='-', scale_div1=1., scale_div2=1., col_posdef=None, ylognorm=False):
    """
    Plot the input distance file consisting of several columns:
    threshold   distance_1  distance_2  ...     distance_M
    ...         ...         ...         ...     ...
    """
    # plot variables
    colors=['b','g']

    # import matrix
    data = np.loadtxt(filein)

    # sort
    idx=np.argsort(data[:,col_thres])
    data=data[idx]

    # set coordinates for plot
    X = data[:,col_thres]
    Y1 = data[:,col_dist1] / scale_div1
    if (col_dist2 != None):
        Y2 = data[:,col_dist2] / scale_div2
    else:
        Y2=None
    if (col_posdef != None):
        Z = data[:,col_posdef]

    # find (last) minimum
    kopt=np.argmin(Y1)
    #kopt=last_minimum(Y)
    xopt=X[kopt]
    yopt=Y1[kopt]
    #print "optimum at xopt={:.2f}  yopt={:.2e}".format(xopt,yopt)

    # create figure
    # plot
    len_xaxis,len_yaxis = 3.5,3.5 #fix here your numbers
    xspace, yspace = .9, .9 # change the size of the void border here.
    x_fig,y_fig = len_xaxis / xspace, len_yaxis / yspace
    figsize=(x_fig,y_fig)
    fig = plt.figure(num='none', facecolor='w',figsize=figsize)
    ax1 = fig.gca()
    if (col_dist2 != None):
        ax2=ax1.twinx()
    else:
        ax2=None

    # plot figure
    ax1.plot(X,Y1,linestyle=linestyle,marker='o', ms=ms, color=colors[0],lw=lw)
    if (annot):
        ax1.plot([xopt],[yopt], 'o', marker='s', ms=2*ms, color=colors[0])
    ax1.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(base=ticks_dx))

    if (ylognorm):
        ax1.set_yscale('log')
    else:
        ax1.yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(base=ticks_dy))

    # labels and axis
    ax1.set_xlabel(xlabel, fontsize='large',labelpad=10)
    if (ax2 == None):
        ax1.set_ylabel(ylabel1,fontsize='large',labelpad=10)
    else:
        ax1.set_ylabel(ylabel1,fontsize='large',color=colors[0],labelpad=10)

        ax2.plot(X,Y2,linestyle=linestyle, marker='v', ms=ms, color=colors[1],lw=lw)
        ax2.set_ylabel(ylabel2,fontsize='large',color=colors[1],labelpad=10)
        if (ylognorm):
            ax2.set_yscale('log')
        else:
            ax2.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(nbins=nticks_max))

    # non-positive definite area
    if (col_posdef != None):
        idx = (Z != 0.)
        xlo = np.min(X[idx])
        xhi = np.max(X[idx])
        ax1.axvspan(xmin=xlo,xmax=xhi,color='red',alpha=0.5, lw=0, label="unphysical")
        ax1.legend(loc='best', fontsize='medium', frameon=False)

    # remove side axis
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.get_xaxis().tick_bottom()
    ax1.get_yaxis().tick_left()
    if (ax2 != None):
        ax2.spines['top'].set_visible(False)
        ax2.spines['left'].set_visible(False)
        ax2.get_xaxis().tick_bottom()
        ax2.get_yaxis().tick_right()

    # make arrow
    if (annot):
        #text="optimal for $\\xi^{{opt}}={:.2f}$".format(xopt)
        text="optimal threshold"
        xy=np.array([xopt,yopt])
        xytext = np.array([0.5,0.9])
        datatr = ax1.transData.inverted()
        axtr = ax1.transAxes
        xytp = axtr.transform(xytext)
        xytext = datatr.transform(xytp)
        if (xy[0] < xytext[0]):
            sgn=+1
        else:
            sgn=-1
        ax1.annotate(text,xy=xy,xycoords='data',xytext=xytext,
                textcoords="data",\
                fontsize='medium',va="center",ha="center",\
                arrowprops=dict(arrowstyle='-|>', connectionstyle='arc3,rad=%.1f' %(sgn*0.2),
                    facecolor='k', shrinkB=5.00))

    # write output and exit
    fileout = os.path.splitext(filein)[0]
    if (ylognorm):
        fileout += "_ylognorm"
    filename=fileout
    for ext in '.png', '.pdf', '.svg':
        fileout = filename + ext
        fig.savefig(fileout,bbox_inches='tight',pad_inches=0,dpi=mydpi)
        print "{:<20s}{:<s}".format("fileout",fileout)
    plt.close('all')
    return fileout

def make_distance_superimposition_plot(distfiles, norms, nij_max, fig_width, fig_height, fontsize, nticks_max, dpi, col_thres, col_dist, xlabel, ylabel, tdir, xbase=None, lw=0.5, ms=2, linestyle='-', scale_div=1.0, ylognorm=False):
    """
    Plot the input distance file consisting of several columns:
    threshold   distance_1  distance_2  ...     distance_M
    ...         ...         ...         ...     ...

    The curves are overlaid.
    """

    # plot variables
    fs_tex=1.15*fontsize
    cmap = matplotlib.cm.viridis

    # sort files
    idx = np.argsort(norms)
    norms=np.array(norms)[idx]
    distfiles = np.array(distfiles)[idx]
    nfiles = len(distfiles)

    # import matrix
    distances = []
    for distfile in distfiles:
        data_full = np.loadtxt(distfile)
        data = []
        data.append(data_full[:,col_thres])
        data.append(data_full[:,col_dist]/scale_div)
        distances.append(np.array(data))

    # create figure
    fig=plt.figure(num='none',facecolor='white', figsize=(fig_width,fig_height))
    ax=fig.gca()
    i0=0
    y0=99.9e99
    # fill the figure
    for i in range(nfiles):
        ## data
        norm=norms[i]
        print "{:<20s}{:.4f}".format("norm", norm)
        distance=distances[i]
        X,Y = distance
        rmax = float(norm) / nij_max

        ## minimum
        kopt=np.argmin(Y)
        xopt=X[kopt]
        yopt=Y[kopt]
        if (yopt < y0):
            y0=yopt
            i0=i

        ## plot curve
        if (nfiles == 1):
            color='b'
        else:
            color=cmap(float(i)/(nfiles-1))

        #ax.plot(X,Y,linestyle=linestyle,marker='o', ms=ms, color=color,lw=lw, label="$\mathcal{{N}} = {:.2g} n_{{max}}$".format(rmax))
        ax.plot(X,Y,linestyle=linestyle,marker='o', ms=ms, color=color,lw=lw, label="$N_c = {:,d}$".format(int(norm)))
        ax.plot([xopt],[yopt], marker='s', ms=2*ms, color=color)

    # print optimum
    print "Minimum for norm={:.1f}".format(norms[i0])
    # labels and axis
    ax.set_xlabel(xlabel, fontsize='large')
    ax.set_ylabel(ylabel,fontsize='large')

    if (xbase == None):
        ax.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(nbins=nticks_max))
    else:
        ax.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(base=xbase))

    if (ylognorm):
        ax.set_yscale('log')
        #ax.set_ylim(0.01,None)
    else:
        ax.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(nbins=nticks_max))


    # legend
    ax.legend(loc='upper center',fontsize='medium', numpoints=3, ncol=3, bbox_to_anchor=(0.5,-0.20), frameon=False)

    # annotate
    ax.annotate("$\max(n_{{ij}})={:,.0f}$".format(np.float_(nij_max)), xy=(1,0), xycoords='axes fraction', va='bottom', ha='right', fontsize='medium')

    # remove side axis
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()

    # write output and exit
    fig.tight_layout(pad=0.1)
    fileout = "distance_superimposition"
    fileout = os.path.join(tdir,fileout)
    if (ylognorm):
        fileout+="_ylognorm"
    filename=fileout

    for ext in '.png', '.pdf', '.svg':
        fileout = filename + ext
        fig.savefig(fileout,bbox_inches='tight',pad_inches=0,dpi=mydpi)
        print "{:<20s}{:<s}".format("fileout",fileout)
    plt.close('all')
    return fileout

