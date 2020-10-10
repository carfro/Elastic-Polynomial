#import os, getopt, sys, math, numpy as np, glob, argparse, re
from numpy import linspace,unique,array,tile,nonzero,meshgrid,argmin,transpose,arange,amax,amin
from math import sqrt,ceil,floor
from scipy import *
from scipy.optimize import leastsq,curve_fit
from lmfit import Model
from ordered_set import OrderedSet
import seaborn as sns
import pandas as pd
from table_generator import exp_val


from matplotlib.widgets import LassoSelector
from matplotlib.path import Path
import matplotlib.pyplot as plt
import matplotlib
from   matplotlib import cm
from   mpl_toolkits.mplot3d import Axes3D
matplotlib.use('TkAgg')
matplotlib.rcParams['text.usetex'] = True

def phi4(x,a,b,c,d,e):
    xu=x
    y=a*xu**4+b*xu**3+c*xu**2+b*xu+e
    return y

def quadfit(x,y):
    cp=[] ; la,lb,lc = polyfit(x,y,2)
    cp.append(la)               
    cp.append(lb) 
    cp.append(lc)
    return cp

def checkfit(x,y): 
    cp=[] ; la,lb,lc,ld = polyfit(x,y,3)
    cp.append(la)               
    cp.append(lb) 
    cp.append(lc)
    cp.append(ld)
    return cp

def sec_order_fit(x,y):
    p=[] ; iia,iib,iic = polyfit(x,y,2)
    p.append(-iib/(2*iia))               	# v0
    p.append(iia*p[0]**2 + iib*p[0] + iic) 	# e0
    p.append((2*iia*p[0]))             		# b0
    p.append(2.)                     		# bp
    return p

def two_var_4order_poly_fit(x,y):
    Pa = polyfit(x[0],y,4)
    Pc = polyfit(x[1],y,4)
    return Pa,Pc

def fit_form_2var_4order(x,k,a0,c0,m11,m22,m12,gam111,gam222,gam112,gam122,\
        del1111,del2222,del1112,del1222,del1122):
    """Two variable polyonomial of deg=2

    Parameters
    ----------
    :x:  TODO
    :k:  TODO
    :a0:
    :c0:
    :m11:
    :m22:
    :m12:
    :gam111:
    :gam222:
    :gam112:
    :gam122:
    :del1111:
    :del2222:
    :del1112:
    :del1222:
    :del1122:

    Returns
    -------
    :returns: Form of two variable fit to fourth order of the energy 

    """
    #return k + 1/2*( m11*(x[0]-a0)**2 + m22*(x[1]-c0)**2 + 2*m12*(x[0]-a0)*(x[1]-c0) ) +\
    #        1/6*( gam111*(x[0]-a0)**3 + gam222*(x[1]-c0)**3 + 3*( gam112*(x[0]-a0)**2*(x[1]-c0) + gam122*(x[0]-a0)*(x[1]-c0)**2 ) ) +\
    #        1/24*( del1111*(x[0]-a0)**4 + del2222*(x[1]-c0)**4 + 4*( del1112*(x[0]-a0)**3*(x[1]-c0)+ del1222*(x[0]-a0)*(x[1]-c0)**3 ) +6*del1122*(x[0]-a0)**2*(x[1]-c0)**2 )
    return k + 1/2*( m11*(x[0]-a0)**2 + m22*(x[1]-c0)**2 + 2*m12*(x[0]-a0)*(x[1]-c0) ) +\
            1/6*( gam111*(x[0]-a0)**3 + gam222*(x[1]-c0)**3 + 3*( gam112*(x[0]-a0)**2*(x[1]-c0) + gam122*(x[0]-a0)*(x[1]-c0)**2 ) ) +\
            1*( del1111*(x[0]-a0)**4 + del2222*(x[1]-c0)**4 +  del1112*(x[0]-a0)**3*(x[1]-c0)+ del1222*(x[0]-a0)*(x[1]-c0)**3 + 1*del1122*(x[0]-a0)**2*(x[1]-c0)**2 )

def E_der_2D(E_poly,x):
    """Calculates the energy derivative of a two dimensional structure.

    :E_poly: TODO
    :param_list: TODO
    :returns: TODO

    """
    m11      = E_poly.params['m11'].value
    m22      = E_poly.params['m22'].value
    m12      = E_poly.params['m12'].value
    gam111   = E_poly.params['gam111'].value
    gam222   = E_poly.params['gam222'].value
    gam112   = E_poly.params['gam112'].value
    gam122   = E_poly.params['gam122'].value
    del1111  = E_poly.params['del1111'].value
    del2222  = E_poly.params['del2222'].value
    del1112  = E_poly.params['del1112'].value
    del1122  = E_poly.params['del1122'].value
    del1222  = E_poly.params['del1222'].value
    a0  = E_poly.params['a0'].value
    c0  = E_poly.params['c0'].value


    E_d1=   1/2*( 2*m11*(x[0]-a0) + 2*m12*(x[1]-c0) ) +\
            1/6*( 3*gam111*(x[0]-a0)**2 + 6*( 2*gam112*(x[0]-a0)*(x[1]-c0) + gam122*(x[1]-c0)**2 ) ) +\
            1*( 4*del1111*(x[0]-a0)**3 +  3*del1112*(x[0]-a0)**2*(x[1]-c0)+ del1222*(x[1]-c0)**3 +2*del1122*(x[0]-a0)*(x[1]-c0)**2 )

    E_d1=-1*a0*E_d1

    E_d2=   1/2*( 2*m22*(x[0]-a0) + 2*m12*(x[0]-c0) ) +\
            1/6*( 3*gam222*(x[1]-a0)**2 + 6*( 2*gam122*(x[0]-a0)*(x[1]-c0) + gam112*(x[0]-c0)**2 ) ) +\
            1*( 4*del2222*(x[1]-a0)**3 + del1112*(x[0]-a0)**3 + del1222*(x[1]-c0)**3 + 2*del1122*(x[0]-a0)**2*(x[1]-c0) )

    E_d2=-1*c0*E_d2
    
    return [E_d1,E_d2]

def two_var_fit(x,y):
    pmodel      = Model(fit_form_2var_4order)
    Pa_i,Pc_i   = two_var_4order_poly_fit(x,y)

    index_emin = argmin(y)
    emin = y[index_emin]
    amin = x[0][index_emin]
    cmin = x[1][index_emin]

    result      = pmodel.fit(y, x=x, k=emin, a0=amin, c0=cmin, m11=Pa_i[2], m22=Pc_i[2], m12=1, gam111=Pa_i[3], gam222=Pc_i[3], gam112=1, gam122=1,\
            del1111=Pa_i[4], del2222=Pc_i[4], del1112=1, del1122=1, del1222=1)
    #result      = pmodel.fit(y, x=x, k=emin, a0=amin, c0=cmin, m11=Pa_i[2], m22=Pc_i[2], m12=0, gam111=Pa_i[3], gam222=Pc_i[3], gam112=0, gam122=0,\
    #        del1111=Pa_i[4], del2222=Pc_i[4], del1112=0, del1122=0, del1222=0)
    return result

def two_var_Refit(twoD_poly,x,y):
    pmodel      = Model(fit_form_2var_4order)
    #Pa_i,Pc_i   = two_var_4order_poly_fit(x,y)

    #index_emin = argmin(y)
    emin = twoD_poly.params["k"]
    amin = twoD_poly.params["a0"]
    cmin = twoD_poly.params["c0"]

    result      = pmodel.fit(y, x=x, k=emin, a0=amin, c0=cmin, m11=twoD_poly.params["m11"], m22=twoD_poly.params["m22"], m12=twoD_poly.params["m12"], gam111=twoD_poly.params["gam111"], gam222=twoD_poly.params["gam222"],\
            gam112=twoD_poly.params["gam112"], gam122=twoD_poly.params['gam122'], del1111=twoD_poly.params["del1111"], del2222=twoD_poly.params["del2222"], \
            del1112=twoD_poly.params["del1112"], del1122=twoD_poly.params["del1122"], del1222=twoD_poly.params["del1222"])
    return result

def fit_form_1var_4order(x,k,a0,m1,gam1,del1):
    """Two variable polyonomial of deg=2

    Parameters
    ----------
    :x:  TODO
    :k:  TODO
    :m11:
    :gam111:
    :del1111:
    :a0:


    Returns
    -------
    :returns: Two variable fit to fourth order of the energy 

    """
    return k + 1/2*m1*(x-a0)**2 + 1/6*gam1*(x-a0)**3 + 1/24*del1*(x-a0)**4

def one_var_fit(x,y):
    pmodel     = Model(fit_form_1var_4order)
    
    index_emin = y.index(min(y))
    
    emin = y[index_emin]
    amin = x[index_emin]

    yfit = [ a-amin for a in y ]
    Pa_i = polyfit(x,yfit,4)

    result      = pmodel.fit(y, x=x, k=emin, a0=amin, m1=Pa_i[2], gam1=Pa_i[3], del1=Pa_i[4])
    return result

def ameos(v,p):
    """ Murnaghan equation of state
    """

    v0 = p[0] ; e0 = p[1] ; b0 = p[2] ; bp = p[3]  
    vv = (v0/v)**(2./3.)
    ff = vv - 1.
    ee = e0 + 9.*b0*v0/16. * ( ff**3*bp + ( 6.-4.*vv )*ff**2 )
    return ee   

def bmeos(v,p):
    """ Birch equation of state
    """

    v0 = p[0] ; e0 = p[1] ; b0 = p[2] ; bp = p[3] 
    bb = 1./(bp-1)
    vv = v0/v    
    ee = e0 + b0/bp*v0/vv*(bb*pow(vv,bp) +1.)-bb*b0*v0
    return ee    

def residuals(p,e,v):
    return e - ameos(v,p)

def bmfit_1D(x,y):
    p = [0,0,0,0] ; p = sec_order_fit(x,y)
    return leastsq(residuals, p, args=(y, x), full_output=1)

def latexify(fig_width=None, fig_height=None, columns=1):
    """Set up matplotlib's RC params for LaTeX plotting.
    Call this before plotting a figure.

    Parameters
    ----------
    fig_width : float, optional, inches
    fig_height : float,  optional, inches
    columns : {1, 2}
    """

    # code adapted from http://www.scipy.org/Cookbook/Matplotlib/LaTeX_Examples

    # Width and max height in inches for IEEE journals taken from
    # computer.org/cms/Computer.org/Journal%20templates/transactions_art_guide.pdf

    assert(columns in [1,2])

    if fig_width is None:
        #fig_width = 3.39 if columns==1 else 6.9 # width in inches
        fig_width = 6.9 if columns==1 else 3.39 # width in inches

    if fig_height is None:
        golden_mean = (sqrt(5)-1.0)/2.0    # Aesthetic ratio
        fig_height = fig_width*golden_mean # height in inches

    MAX_HEIGHT_INCHES = 8.0
    if fig_height > MAX_HEIGHT_INCHES:
        print("WARNING: fig_height too large:" + fig_height +
              "so will reduce to" + MAX_HEIGHT_INCHES + "inches.")
        fig_height = MAX_HEIGHT_INCHES

    params = {'backend': 'ps',
             # 'text.latex.preamble': ['\usepackage{gensymb}'],
              'axes.labelsize': 12, # fontsize for x and y labels (was 11)
              'axes.titlesize': 12,
              #'text.fontsize': 8, # was 10
              'legend.fontsize': 8, # was 10
              'xtick.labelsize': 12,
              'ytick.labelsize': 12,
              'text.usetex': True,
              'figure.figsize': [fig_width,fig_height],
              'font.family': 'serif'
    }

    matplotlib.rcParams.update(params)
    return

def format_axes(ax):

    for spine in ['top', 'right']:
        ax.spines[spine].set_visible(False)

    for spine in ['left', 'bottom']:
        ax.spines[spine].set_color('k')
        ax.spines[spine].set_linewidth(0.5)

    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')

    for axis in [ax.xaxis, ax.yaxis]:
        axis.set_tick_params(direction='out', color='k')

    return ax

class SelectFromCollection(object):
    """Select indices from a matplotlib collection using `LassoSelector`.

    Selected indices are saved in the `ind` attribute. This tool fades out the
    points that are not part of the selection (i.e., reduces their alpha
    values). If your collection has alpha < 1, this tool will permanently
    alter the alpha values.

    Note that this tool selects collection objects based on their *origins*
    (i.e., `offsets`).

    Parameters
    ----------
    ax : :class:`~matplotlib.axes.Axes`
        Axes to interact with.

    collection : :class:`matplotlib.collections.Collection` subclass
        Collection you want to select from.

    alpha_other : 0 <= float <= 1
        To highlight a selection, this tool sets all selected points to an
        alpha value of 1 and non-selected points to `alpha_other`.
    """

    def __init__(self, ax, collection, alpha_other=0.3):
        self.canvas = ax.figure.canvas
        self.collection = collection
        self.alpha_other = alpha_other

        self.xys = collection.get_offsets()
        self.Npts = len(self.xys)

        # Ensure that we have separate colors for each object
        self.fc = array(collection.get_facecolors())
        if len(self.fc) == 0:
            raise ValueError('Collection must have a facecolor')
        elif len(self.fc) == 1:
            self.fc = tile(self.fc, (self.Npts, 1))

        self.lasso = LassoSelector(ax, onselect=self.onselect)
        self.ind = []

    def onselect(self, verts):
        path = Path(verts)
        self.ind = nonzero(path.contains_points(self.xys))[0]
        self.fc[:, -1] = self.alpha_other
        self.fc[self.ind, -1] = 1
        self.collection.set_facecolors(self.fc)
        self.canvas.draw_idle()

    def disconnect(self):
        self.lasso.disconnect_events()
        self.fc[:, -1] = 1
        self.collection.set_facecolors(self.fc)
        self.canvas.draw_idle()

def data_plot_oneD(result_res_dict, *fmt_list):
    """
    A helper function to make a graph 

    Parameters
    ----------
    param_list : Axes
        The axes to draw to

    v_list : array
       The x data

    e_list : array
       The y data

    *fmt_list : dict
       Dictionary of kwargs to pass to ax.plot

    Returns
    -------
    out : list
        list of artists added
    """
#    plt.
#    out = ax.plot(data1, data2, **param_dict)

    fig_width=10
    fig_height=15

    matplotlib.rcParams.update(matplotlib.rcParamsDefault)

    # FOR LATEX PLOTTING
    matplotlib.rcParams.update({'text.usetex': True,
                      'axes.labelsize': 14, # fontsize for x and y labels (was 11)
                      'axes.titlesize': 14,
                      #'text.fontsize': 8, # was 10
                      'legend.fontsize': 14, # was 10
                      'xtick.labelsize': 14,
                      'ytick.labelsize': 14,
                      'figure.figsize': [fig_width,fig_height],
                      'font.family': 'serif'})

    fig,axs = plt.subplots(3,figsize=(fig_width, fig_height))

    colors = iter(cm.rainbow(linspace(0, 1, len(result_res_dict))))
    min_e=[]
    max_e=[]
    min_p=[]

    for res_dict in result_res_dict:
        v_list=res_dict["v_list"]; e_list=res_dict["e_list"]; param_list=res_dict["param_list"]; lm_poly=res_dict["poly"]; BM_poly=res_dict["BM_poly"]; geo=res_dict["geo"]

        #print(param_list)
        #print(type(res_dict[9]))

        index_emin = e_list.index(min(e_list))
        min_e.append(min(e_list))
        max_e.append(max(e_list))
        min_p.append(param_list[index_emin])

        if (geo!="3L-c"):
            v_grid=linspace(v_list[0],v_list[-1],1000)
        a_grid=linspace(param_list[0],param_list[-1],1000)
        col = next(colors)

        # Plot of data
        axs[0].plot(param_list,e_list,'o',color=col,label=res_dict["name"]+" -"+res_dict["functional"])
        if len(res_dict)>10:
            axs[0].plot([res_dict["exp"],res_dict["exp"]],[min(e_list),max(e_list)],'--k')
            axs[1].plot([res_dict["exp"],res_dict["exp"]],[min(e_list),max(e_list)],'--k')
            axs[2].plot([res_dict["exp"],res_dict["exp"]],[min(e_list),max(e_list)],'--k')

        # Plot of 4th deg poly
        axs[1].plot(param_list,e_list,'o', color=col,label="data")
        #axs[1].plot(param_list ,lm_poly.init_fit, '--',color=col, label="Initial fit")
        axs[1].plot(a_grid,lm_poly.eval(x=a_grid), '-',color=col, label="Best fit")
        # Plot of Murnaghan fit
        if (geo!="3L-c"):
            axs[2].plot(v_list,e_list,'o',color=col,label="data")
            axs[2].plot(v_grid,ameos(v_grid,BM_poly[0]),'-',color=col,label="Murnaghan fit")

    # Plot of data
    axs[0].plot(min_p,min_e,'ko',label="Min values")
    axs[0].set_xlabel(r"a [Å]")
    axs[0].set_ylabel(r"Energy [eV]")
    axs[0].legend(loc="best")
    axs[0].title.set_text('Data plot')
    #axs[0].set_ylim(bottom=(floor(min(min_e)*100)/100.0),top=(ceil(max(max_e)*100)/100.0))

    # Plot of 4th deg poly
    axs[1].set_xlabel(r"a [Å]")
    axs[1].set_ylabel(r"Energy [eV]")
    #axs[1].legend(loc="best")
    axs[1].title.set_text('Fourth order poly. fit')
    #axs[1].set_ylim(bottom=(floor(min(min_e)*100)/100.0),top=(ceil(max(max_e)*100)/100.0))

    # Plot of Murnaghan fit
    if (geo!="3L-c"):
        axs[2].set_xlabel(r"V  [Å$\displaystyle ^3$]")
        axs[2].set_ylabel(r"Energy [eV]")
        axs[2].title.set_text('Murnaghan (volume) fit')
    #axs[2].set_ylim(bottom=(floor(min(min_e)*100)/100.0),top=(ceil(max(max_e)*100)/100.0))

    return fig

#def pdf_plot_oneD(arg_plot,param_list,v_list,e_list,BM_poly,lm_poly, *fmt_list):
def pdf_plot_oneD(arg_plot,result_dict_list, *fmt_list):
    #TODO: FIX SO THAT YOU CAN HAVE .PDF IN NAME IF YOU'D LIKE
    """
    A helper function to make a graph 

    Parameters
    ----------
    param_list : Axes
        The axes to draw to

    v_list : array
       The x data

    e_list : array
       The y data

    *fmt_list : dict
       Dictionary of kwargs to pass to ax.plot

    Returns
    -------
    out : list
        list of artists added
    """
    latexify()

    species_list        = [ res_dict['name'] for res_dict in result_dict_list ]
    print("species_list",species_list)
    functional_list     = [ res_dict['functional'] for res_dict in result_dict_list ]
    print("functional_list",functional_list)
    unique_species      = OrderedSet(species_list)
    print("unique_species",unique_species)
    #unique_functionals  = unique(functional_list).tolist()
    unique_functionals  = OrderedSet(functional_list)
    #unique_functionals.sort(key  = len)
    print("unique_functionals",unique_functionals)

    len_func  = len(unique_functionals)
    len_spec  = len(unique_species)
    size_grid = len_func*len_spec

    index_species       = [[i for i, e in enumerate(species_list) if e == u_s] for u_s in unique_species]
    print("index_species",index_species)
    index_functionals   = [[i for i, e in enumerate(functional_list) if e == u_s] for u_s in unique_functionals]
    print("index_functionals",index_functionals)
    plot_order          = [[ i for j in range(len_func) for i in index_functionals[j] if i in S] for S in index_species]
    print("PLOT ORDER",plot_order)
    

    colors = iter(cm.rainbow(linspace(0, 1, len_spec)))

    if(arg_plot!="NO PDF"):
        name = '%s_poly.pdf'%(arg_plot)
    #fig, ax = plt.subplots(nrows=len_func,ncols=len_spec,figsize=(20,25))
    fig, ax = plt.subplots(nrows=len_func,ncols=len_spec,figsize=(5,10))
    #fig, ax = plt.subplots(nrows=len_func,ncols=len_spec,figsize=(10,15))

    if(len(result_dict_list)>1):
        fig.text(0.5, 0.01, r"{\Large a [Å]}", ha='center', va='center')
        i=0
        j=0
        for plot_list in plot_order:
            col = next(colors)
            for cur_plot in plot_list:
                #print("cur_plot",cur_plot)
                #print("index",i)
                res_dict = result_dict_list[cur_plot]


                v_list=res_dict['v_list']; e_list=res_dict['e_list']; param_list=res_dict['param_list']; lm_poly=res_dict['poly'];BM_poly=res_dict["BM_poly"]

                a_grid=linspace(param_list[0],param_list[-1],1000)
                if ( len_spec==1 ) or ( len_func==1 ):

                    ax[i].plot(param_list,e_list,'o',color=col,label="data")
                    ax[i].plot(a_grid,lm_poly.eval(x=a_grid),'-',color=col, label="Best fit")
                    if len(res_dict)>10:
                        ax[i].axvline(res_dict["exp"],color='k',ls='--',lw=0.5)
                    locs=list(ax[i].get_xticks())
                    locs=[ round(e,3) for e in locs ]
                    labels = ["$"+str(w)+"$" for w in locs]
                    locsy=list(ax[i].get_yticks())
                    locs+=[res_dict["exp"]]
                    labels+=["$"+str(res_dict["exp"])+"$"]
                    ax[i].set_xticks(locs)
                    ax[i].set_xticklabels(labels)
                    ax[i].get_xticklabels()[-1].set_rotation(45)
                    ax[i].get_xticklabels()[-1].set_position((1,-0.010))
                    ax[i].set_yticks(locsy)
                    plt.tight_layout()
                    format_axes(ax[i])
                    ax[i].title.set_text('%s using %s'%(res_dict["name"],res_dict['functional']))
                    #ax[i,j].set_xlabel(r"a [Å]")
                    #ax[i,j].set_ylabel(r"Energy [eV]")
                    ax[i].grid(which='major',linestyle='--',dashes=(5,8), linewidth='0.25')
                    if (i==0) and (j==0):
                        ax[i].set_ylabel(r"{\Large Energy [eV]}")
                    #if (i==ceil(len_spec/2)-1) and (j==0):
                    #    ax[i].set_ylabel(r"{\Large Energy [eV]}")
                    plt.draw()
                else:
                    #print( "index_functionals[i]",index_functionals[i] )
                    #print( "cur_plot",cur_plot)
                    #print( "i",i )
                    if not ( cur_plot in index_functionals[i]):
                       #print()
                       #print("Searching for correct column...")
                       while not ( cur_plot in index_functionals[i]):  
                           if(i >=len(index_functionals)-1):
                                i=0
                                #print("starting search from beginning")
                           #print( "index_functionals[i]",index_functionals[i] )
                           #print( "i",i )
                           i+=1
                    ax[i,j].plot(param_list,e_list,'o',color=col,label="data")
                    ax[i,j].plot(a_grid,lm_poly.eval(x=a_grid),'-',color=col, label="Best fit")
                    if(res_dict["name"]=="MgO"):
                        print("MgO len of e_list",len(e_list))
                    if len(res_dict)>10:
                        ax[i,j].axvline(res_dict["exp"],color='k',ls='--',lw=0.5)
                    locs=list(ax[i,j].get_xticks())
                    locs=[ round(e,3) for e in locs ]
                    labels = ["$"+str(w)+"$" for w in locs]
                    locsy=list(ax[i,j].get_yticks())
                    locs+=[res_dict["exp"]]
                    labels+=["$"+str(res_dict["exp"])+"$"]
                    ax[i,j].set_xticks(locs)
                    ax[i,j].set_xticklabels(labels)
                    ax[i,j].get_xticklabels()[-1].set_rotation(45)
                    ax[i,j].get_xticklabels()[-1].set_position((1,-0.015))
                    ax[i,j].set_yticks(locsy)
                    plt.tight_layout()
                    format_axes(ax[i,j])
                    ax[i,j].title.set_text('%s using %s'%(res_dict["name"],res_dict['functional']))
                    #ax[i,j].set_xlabel(r"a [Å]")
                    #ax[i,j].set_ylabel(r"Energy [eV]")
                    ax[i,j].grid(which='major',linestyle='--',dashes=(5,8), linewidth='0.25')
                    if (i==ceil(len_func/2)-1) and (j==0):
                        ax[i,j].set_ylabel(r"{\Large Energy [eV]}")
                    plt.draw()
                i+=1
                #print("END OF LOOOP")
            j+=1
            i = i if ( len_spec==1 ) or ( len_func==1 ) else 0 # width in inches
    else:
        res_dict = result_dict_list[0]
        v_list=res_dict['v_list']; e_list=res_dict['e_list']; param_list=res_dict['param_list']; lm_poly=res_dict['poly'];BM_poly=res_dict["BM_poly"]

        ax.grid(which='major',linestyle='--',dashes=(5,8), linewidth='0.25')
        ax.set_ylabel(r"{Energy [eV]}")
        ax.set_xlabel(r"{c [Å]}")

        a_grid=linspace(param_list[0],param_list[-1],1000)
        ax.plot(param_list,e_list,'o',color="blue",label="data")
        ax.plot(a_grid,lm_poly.eval(x=a_grid),'-',color="blue", label="Best fit")

        locs=list(ax.get_xticks())
        locs=[ round(e,3) for e in locs ]
        labels = ["$"+str(w)+"$" for w in locs]
        locsy=list(ax.get_yticks())

        ax.set_xticks(locs)
        ax.set_xticklabels(labels)
        #ax.get_xticklabels()[-1].set_rotation(45)
        #ax.get_xticklabels()[-1].set_position((1,-0.010))
        ax.set_yticks(locsy)
        plt.tight_layout()
        format_axes(ax)


        plt.draw()
        plt.show()

    if(arg_plot!="NO PDF"):
        fig.savefig(name,format="pdf")
        plt.show()
        plt.close(name)
        print('Plot ',name,' was created.\n')
        return
    else:
        return fig

def set_violin_prop(violin,i):
#   for partname in ('cmins','cmaxes'):
#        vp = violin[partname]
#        vp.set_edgecolor(colors[i])
#        vp.set_linewidth(1)
#        vp.set_color(colors[i])
    colors = array(['C0', 'C1', 'C3', 'C2', 'C4', 'C5', 'k'])
    violin['cmedians'].set_color('black')
    violin['cmedians'].set_linewidth(2)
    for parts in violin['bodies']:
        parts.set_color(colors[i])

def pdf_violin_plot_oneD(arg_plot,result_dict_list, *fmt_list):
    #TODO: FIX SO THAT YOU CAN HAVE .PDF IN NAME IF YOU'D LIKE
    """
    A helper function to make a graph 

    Parameters
    ----------
    param_list : Axes
        The axes to draw to

    v_list : array
       The x data

    e_list : array
       The y data

    *fmt_list : dict
       Dictionary of kwargs to pass to ax.plot

    Returns
    -------
    out : list
        list of artists added
    """
    latexify()

    fig_width=10
    fig_height=15

    matplotlib.rcParams.update(matplotlib.rcParamsDefault)

    # FOR LATEX PLOTTING
    matplotlib.rcParams.update({'text.usetex': True,
                      'axes.labelsize': 28, # fontsize for x and y labels (was 11)
                      'axes.titlesize': 28,
                      #'text.fontsize': 8, # was 10
                      'legend.fontsize': 28, # was 10
                      'xtick.labelsize': 28,
                      'ytick.labelsize': 28,
                      'figure.figsize': [fig_width,fig_height],
                      'font.family': 'serif'})

    species_list        = [ res_dict['name'] for res_dict in result_dict_list ]
    #print("species_list",species_list)

    functional_list     = [ res_dict['functional'] for res_dict in result_dict_list ]
    #print("functional_list",functional_list)

    unique_species      = OrderedSet(species_list)

    unique_functionals = OrderedSet(functional_list)
    print("ORDERED SET",unique_functionals)
    #print()

    len_func  = len(unique_functionals)

    #print()
    a_0=[[res_dict["res_poly"][0] for res_dict in result_dict_list if res_dict['functional'] == u_s] for u_s in unique_functionals]
    #print(a_0)
    for a_sub in a_0:
        j=0
        for spec in unique_species:
            a_sub[j] = (a_sub[j]-exp_val(spec,0))
            #a_sub[j] = (a_sub[j])
            j+=1
    #print()
    #print(a_0)
    #print()

    #print()
    e_coh=[[res_dict["res_poly"][1] for res_dict in result_dict_list if res_dict['functional'] == u_s] for u_s in unique_functionals] 
    #print("E_coh",e_coh)
    for e_sub in e_coh:
        j=0
        for spec in unique_species:
            #print("system",spec)
            #print("calc value",e_sub[j])
            #print("exp value",exp_val(spec,2))
            e_sub[j] = (e_sub[j]-exp_val(spec,1))/(exp_val(spec,1))*100
            j+=1
    #print()
    #print("(E_coh-E_coh)/E_coh",e_coh)
    #print()
    #exit()

    #print()
    b_0=[[res_dict["res_poly"][2] for res_dict in result_dict_list if res_dict['functional'] == u_s] for u_s in unique_functionals] 
    #print("B_0",b_0)
    for b_sub in b_0:
        j=0
        for spec in unique_species:
            #print("system",spec)
            #print("calc value",b_sub[j])
            #print("exp value",exp_val(spec,2))
            b_sub[j] = (b_sub[j]-exp_val(spec,2))/(exp_val(spec,2))*100
            j+=1
    #print()
    #print("(B_0-B_ref)/B_ref",b_0)
    #print()
    #exit()

    colors = iter(cm.rainbow(linspace(0, 1, len_func)))

    name = '%s_VIOLIN.pdf'%(arg_plot)
    fig, axs = plt.subplots(nrows=2,ncols=1,figsize=(fig_width,fig_height))

    def plot_data(start,result,ax,labels=[],fact=1):
        for i in range(len_func):
            bp = ax.boxplot(result[i]*fact,positions=[start+i], widths=[0.2], patch_artist = True,
                    showmeans=True,
                    showfliers=True)

            outliers = [flier.get_ydata() for flier in bp["fliers"]]
            #for iX40, val in enumerate(result[i]*fact):
            #    if len(outliers):
            #       for out in outliers[0]:
            #           if abs(out - val) < 1e-10:
            #                #print(unique_functionals[i],iX40, val)
            plt.setp(bp['medians'], linewidth=2,color='k',alpha=1.0)
            plt.setp(bp['boxes'], linewidth=1,color='k',alpha=0.5)
            plt.setp(bp['fliers'], marker='o', markeredgecolor='k', markersize=5, markerfacecolor='white', markeredgewidth=1.0, lw=1,alpha=0.7)
            plt.setp(bp['means'], marker='D',markersize=4, markerfacecolor='white', markeredgecolor='white',alpha=1.0)
            violin = ax.violinplot(result[i]*fact,positions=[start+i], widths=[0.8],
                               showextrema=False, showmeans=False, showmedians=True)

            set_violin_prop(violin,i)
        ax.set_xlim(-0.5,3.5)
        ax.axhline(ls = '--', color='gray',lw=2)

    plot_data(0,a_0,axs[0],fact=100)
    plot_data(0,e_coh,axs[1],fact=100)
    #plot_data(0,b_0,axs[2],fact=100)

    axs[0].set_ylabel(r'$ \Delta a \, ({\rm \AA})$')
    axs[1].set_ylabel(r'$ (E- E_{\rm ref})/ E_{\rm ref}\,(\%)$')
    #axs[2].set_ylabel(r'$ (B_0- B_{0,\, {\rm ref}})/ B_{0,\, {\rm ref}}\,(\%)$')

    axs[0].set_xticks(list(range(len_func)))
    axs[0].set_xticklabels([])

    axs[1].set_xticks(list(range(len_func)))
    #axs[1].set_xticklabels([])
    axs[1].set_xticklabels(unique_functionals,rotation=-60)

    #axs[2].set_xticks(list(range(len_func)))
    #axs[2].set_xticklabels(unique_functionals,rotation=-60)
    
    plt.tight_layout()

    plt.show()
    fig.savefig(name,format="pdf")
    plt.close()
    print('Plot ',name,' was created.\n')

    return 

def pdf_bar_plot_oneD(arg_plot,result_dict_list, *fmt_list):
    #TODO: FIX SO THAT YOU CAN HAVE .PDF IN NAME IF YOU'D LIKE
    """
    A helper function to make a graph 

    Parameters
    ----------
    param_list : Axes
        The axes to draw to

    v_list : array
       The x data

    e_list : array
       The y data

    *fmt_list : dict
       Dictionary of kwargs to pass to ax.plot

    Returns
    -------
    out : list
        list of artists added
    """
    latexify()

    fig_width=15
    fig_height=20

    matplotlib.rcParams.update(matplotlib.rcParamsDefault)

    # FOR LATEX PLOTTING
    matplotlib.rcParams.update({'text.usetex': True,
                      'axes.labelsize': 18, # fontsize for x and y labels (was 11)
                      'axes.titlesize': 18,
                      #'text.fontsize': 8, # was 10
                      'legend.fontsize': 18, # was 10
                      'xtick.labelsize': 18,
                      'ytick.labelsize': 18,
                      'figure.figsize': [fig_width,fig_height],
                      'font.family': 'serif'})

    species_list        = [ res_dict['name'] for res_dict in result_dict_list ]
    print("species_list",species_list)

    functional_list     = [ res_dict['functional'] for res_dict in result_dict_list ]
    print("functional_list",functional_list)

    unique_species      = OrderedSet(species_list)
    print("unique_species",unique_species)

    unique_functionals = OrderedSet(functional_list)
    print("unique_functionals",unique_functionals)
    #print()

    len_func  = len(unique_functionals)
    
    #print()
    a_0=[[res_dict["res_poly"][0] for res_dict in result_dict_list if res_dict['functional'] == u_s] for u_s in unique_functionals]
    #print(a_0)
    for a_sub in a_0:
        j=0
        for spec in unique_species:
            a_sub[j] = (a_sub[j]-exp_val(spec,0))/(exp_val(spec,0))*100
            j+=1
    #print()
    #print(a_0)

    e_coh=[[res_dict["res_poly"][1] for res_dict in result_dict_list if res_dict['functional'] == u_s] for u_s in unique_functionals] 
    for e_sub in e_coh:
        j=0
        for spec in unique_species:
            e_sub[j] = (e_sub[j]-exp_val(spec,1))/(exp_val(spec,1))*100
            j+=1

    b_0=[[res_dict["res_poly"][2] for res_dict in result_dict_list if res_dict['functional'] == u_s] for u_s in unique_functionals] 
    for b_sub in b_0:
        j=0
        for spec in unique_species:
            b_sub[j] = (b_sub[j]-exp_val(spec,2))/(exp_val(spec,2))*100
            j+=1

    func_MAD=[]
    for i, func in enumerate(unique_functionals):
        func_MAD.append( [mean(a_0[i]),mean(e_coh[i]),mean(b_0[i])] )

    #colors = iter(cm.rainbow(linspace(0, 1, len_func)))
    colors = array(['C0', 'C1', 'C3', 'C2', 'C4', 'C5', 'k'])

    barWidth = 0.40
    midFac = (len(func_MAD[0])/2-0.5)
    name = '%s_BAR.pdf'%(arg_plot)
    rows=math.floor(len_func/2.)
    cols=math.ceil(len_func/2.)
    fig, axs = plt.subplots(nrows=rows,ncols=cols,figsize=(fig_width,fig_height),sharex='col', sharey='row',gridspec_kw={'hspace': 0, 'wspace': 0})

    ##fig.text(0.5, 0.01, r"{\Large a [Å]}", ha='center', va='center')

    ind_MAD=0
    labels=['$a_0$','$E_{coh}$','$B_0$']
    for i in range(cols):
        for j in range(rows):
            print("row",2-j,"col",i+1)
            print("ind MAD",ind_MAD)
            print(unique_functionals[ind_MAD])
            r = arange(len(func_MAD[ind_MAD]))
            print(func_MAD[ind_MAD])
            col = colors[ind_MAD]
            axs[1-j][i].grid(which='major',linestyle='--',dashes=(5,8), linewidth='0.4')
            axs[1-j][i].bar(r, func_MAD[ind_MAD], color=col, width=barWidth, edgecolor="white")
            axs[1-j][i].set_ylim([2*round(amin(func_MAD)/2),2*round(amax(func_MAD)/2)])
            axs[1-j][i].axhline(ls = '--', color='gray',lw=1.5)

            if (i==0) and (j==1):
                locs_y=list(axs[j][i].get_yticks())
                axs[1-j][i].set_yticks([e for e in locs_y[1:]])
            #elif :
            #    axs[i][j].bar(r, func_MAD[ind_MAD], color=col, width=barWidth, edgecolor="white",sharex = prev_axs)
            #    prev_axs=axs[i][j]
            r=[x + barWidth for x in r]
            if (j!=rows-1):
                axs[1-j][i].set_xticks([r + midFac*barWidth for r in range(len(func_MAD[ind_MAD]))])
                axs[1-j][i].set_xticklabels([])
            else:
                axs[1-j][i].set_xticks([r for r in range(len(func_MAD[ind_MAD]))])
                axs[1-j][i].set_xticklabels(labels)
            ind_MAD+=1
            #plt.legend(loc='upper left')
    fig.text(0.04, 0.5, 'Mean relative deviation $(\%)$', va='center', rotation='vertical',size=24)
    #plt.setp(axs, ylim=[amin(func_MAD),amax(func_MAD)],yticks=[amin(func_MAD),amax(func_MAD)])

    plt.show()
    fig.savefig(name,format="pdf")
    plt.close(name)
    print('Plot ',name,' was created.\n')

    return 

def pdf_violin_plot_oneD_OLD(arg_plot,result_dict_list, *fmt_list):
    #TODO: FIX SO THAT YOU CAN HAVE .PDF IN NAME IF YOU'D LIKE
    """
    A helper function to make a graph 

    Parameters
    ----------
    param_list : Axes
        The axes to draw to

    v_list : array
       The x data

    e_list : array
       The y data

    *fmt_list : dict
       Dictionary of kwargs to pass to ax.plot

    Returns
    -------
    out : list
        list of artists added
    """
    latexify()

    fig_width=15
    fig_height=20

    matplotlib.rcParams.update(matplotlib.rcParamsDefault)

    # FOR LATEX PLOTTING
    matplotlib.rcParams.update({'text.usetex': True,
                      'axes.labelsize': 16, # fontsize for x and y labels (was 11)
                      'axes.titlesize': 16,
                      #'text.fontsize': 8, # was 10
                      'legend.fontsize': 16, # was 10
                      'xtick.labelsize': 16,
                      'ytick.labelsize': 16,
                      'figure.figsize': [fig_width,fig_height],
                      'font.family': 'serif'})

    species_list        = [ res_dict['name'] for res_dict in result_dict_list ]
    #print("species_list",species_list)

    functional_list     = [ res_dict['functional'] for res_dict in result_dict_list ]
    #print("functional_list",functional_list)

    unique_species      = unique(species_list)
    #print("unique_species",unique_species)

    unique_functionals = OrderedSet(functional_list)
    print("ORDERED SET",unique_functionals)
    #print()

    len_func  = len(unique_functionals)
    
    #ene=[[res_dict["res_poly"][0] for res_dict in result_dict_list if res_dict['functional'] == u_s] for u_s in unique_functionals] 
    a_0=[[res_dict["res_poly"][0] for res_dict in result_dict_list if res_dict['functional'] == u_s] for u_s in unique_functionals]
    #a_tuples=[ (a_0[i][j],f) for i in enum]
    e_coh=[[res_dict["res_poly"][1] for res_dict in result_dict_list if res_dict['functional'] == u_s] for u_s in unique_functionals] 
    b_0=[[res_dict["res_poly"][2] for res_dict in result_dict_list if res_dict['functional'] == u_s] for u_s in unique_functionals] 
    #print(a_0)
    #print(e_coh)
    #print(b_0)

    a_tuple_list=[]
    i=0
    for a_sub in a_0:
        a_tuple_list.extend((list(zip(a_sub,[unique_functionals[i] for j in enumerate(a_sub)]))))
        i+=1
    print(a_tuple_list)
    a_data=pd.DataFrame(a_tuple_list,columns=["$a_0$","functional"])
    print(a_data)
        
    e_tuple_list=[]
    i=0
    for e_sub in e_coh:
        e_tuple_list.extend((list(zip(e_sub,[unique_functionals[i] for j in enumerate(e_sub)]))))
        i+=1
    print(e_tuple_list)
    e_data=pd.DataFrame(e_tuple_list,columns=["$E_{coh}$","functional"])
    print(e_data)

    b_tuple_list=[]
    i=0
    for b_sub in b_0:
        b_tuple_list.extend((list(zip(b_sub,[unique_functionals[i] for j in enumerate(b_sub)]))))
        i+=1
    print(b_tuple_list)
    b_data=pd.DataFrame(b_tuple_list,columns=["$B_0$","functional"])
    print(b_data)

    #atuple_data=pd.DataFrame(a_tuple,columns=["values","functional"])
    #print(atuple_data)
    #print()
    #
    #a_0=transpose(array(a_0))
    #a_data = pd.DataFrame(a_0,columns=unique_functionals)
    #print(a_data)
    ##print(data)
    #print()
    #tips = sns.load_dataset("tips")
    #print(tips)
    #print()

    colors = iter(cm.rainbow(linspace(0, 1, len_func)))

    #name = 'TEST_VIOLIN.pdf'
    name = '%s_VIOLIN.pdf'%(arg_plot)
    fig, axs = plt.subplots(nrows=3,ncols=1,figsize=(fig_width,fig_height))
    #fig1, axs1 = plt.subplots(nrows=3,ncols=1,figsize=(20,25))
    #fig2, axs2 = plt.subplots(nrows=3,ncols=1,figsize=(20,25))

    #fig.text(0.5, 0.01, r"{\Large a [Å]}", ha='center', va='center')

    #bp0 = axs[0].violinplot(a_0)
    #sns.set(style="whitegrid")
    sns.violinplot(x="functional",y="$a_0$",data=a_data,ax=axs[0])
    #sns.violinplot(x="functional",y="$a_0$",data=a_data,ax=axs[0], inner="stick", bw=.2, cut=0)
    sns.violinplot(x="functional",y="$a_0$",data=a_data,ax=axs[1], bw=.2, inner="quartile")
    #axs[2].violinplot(a_0)
    sns.violinplot(x="functional",y="$a_0$",data=a_data,ax=axs[2], bw=.2, scale="count")
    #axs[0].set_ylabel(r"$a$ [Å]")
    #axs[0].set_xlabel(r"Functional")
    #axs[0].title.set_text('a_0')


    #bp1 = axs[1].violinplot(e_coh)
    #sns.violinplot(x="functional",y="$E_{coh}$",data=e_data,ax=axs[1], inner="stick",bw=.2, scale="width")
    #axs[1].set_ylabel(r"$E_{coh}$ [eV]")
    #axs[1].set_xlabel(r"Functional")
    ##axs[1].title.set_text('E_coh')


    #bp2 = axs[2].violinplot(b_0)
    #sns.violinplot(x="functional",y="$B_0$",data=b_data,ax=axs[2], inner="stick",bw=.2, scale="count")
    #axs[2].set_ylabel(r"$B_0$ [GPa]")
    #axs[2].set_xlabel(r"Functional")
    #axs[2].title.set_text('B_0')


    plt.show()
    fig.savefig(name,format="pdf")
    plt.close(name)
    print('Plot ',name,' was created.\n')

    return 

def pdf_bar_plot_oneD_OLD(arg_plot,result_dict_list, *fmt_list):
    #TODO: FIX SO THAT YOU CAN HAVE .PDF IN NAME IF YOU'D LIKE
    """
    A helper function to make a graph 

    Parameters
    ----------
    param_list : Axes
        The axes to draw to

    v_list : array
       The x data

    e_list : array
       The y data

    *fmt_list : dict
       Dictionary of kwargs to pass to ax.plot

    Returns
    -------
    out : list
        list of artists added
    """
    latexify()

    fig_width=12
    fig_height=15

    matplotlib.rcParams.update(matplotlib.rcParamsDefault)

    # FOR LATEX PLOTTING
    matplotlib.rcParams.update({'text.usetex': True,
                      'axes.labelsize': 20, # fontsize for x and y labels (was 11)
                      'axes.titlesize': 22,
                      #'text.fontsize': 8, # was 10
                      'legend.fontsize': 20, # was 10
                      'xtick.labelsize': 20,
                      'ytick.labelsize': 20,
                      'figure.figsize': [fig_width,fig_height],
                      'font.family': 'serif'})

    species_list        = [ res_dict['name'] for res_dict in result_dict_list ]
    print("species_list",species_list)

    functional_list     = [ res_dict['functional'] for res_dict in result_dict_list ]
    print("functional_list",functional_list)

    unique_species      = OrderedSet(species_list)
    print("unique_species",unique_species)

    unique_functionals = OrderedSet(functional_list)
    print("unique_functionals",unique_functionals)
    #print()

    len_func  = len(unique_functionals)
    
    #print()
    a_0=[[res_dict["res_poly"][0] for res_dict in result_dict_list if res_dict['functional'] == u_s] for u_s in unique_functionals]
    #print(a_0)
    for a_sub in a_0:
        j=0
        for spec in unique_species:
            a_sub[j] = (a_sub[j]-exp_val(spec,0))/(exp_val(spec,0))*100
            j+=1
    #print()
    #print(a_0)

    e_coh=[[res_dict["res_poly"][1] for res_dict in result_dict_list if res_dict['functional'] == u_s] for u_s in unique_functionals] 
    for e_sub in e_coh:
        j=0
        for spec in unique_species:
            e_sub[j] = (e_sub[j]-exp_val(spec,1))/(exp_val(spec,1))*100
            j+=1

    b_0=[[res_dict["res_poly"][2] for res_dict in result_dict_list if res_dict['functional'] == u_s] for u_s in unique_functionals] 
    for b_sub in b_0:
        j=0
        for spec in unique_species:
            b_sub[j] = (b_sub[j]-exp_val(spec,2))/(exp_val(spec,2))*100
            j+=1

    colors = array(['C0', 'C1', 'C3', 'C2', 'C4', 'C5', 'k'])

    barWidth = 0.15
    midFac = (len(unique_functionals)/2-0.5)
    name = '%s_BAR.pdf'%(arg_plot)
    fig, axs = plt.subplots(nrows=2,ncols=1,figsize=(fig_width,fig_height))

    ##fig.text(0.5, 0.01, r"{\Large a [Å]}", ha='center', va='center')

    r = arange(len(a_0[0]))
    i=0
    for a_sub in a_0:
        #col = next(colors)
        col = colors[i]
        axs[0].bar(r, a_sub, color=col, width=barWidth, edgecolor="white", label=unique_functionals[i])
        r=[x + barWidth for x in r]
        i+=1
    axs[0].set_xticks([r + midFac*barWidth for r in range(len(a_0[0]))])
    axs[0].set_xticklabels(unique_species)
    axs[0].set_ylabel("Rel. error in $a_0$\, (\%)")
    #plt.legend(loc='upper left')


    r = arange(len(e_coh[0]))
    i=0
    for e_sub in e_coh:
        #col = next(colors)
        col = colors[i]
        axs[1].bar(r, e_sub, color=col, width=barWidth, edgecolor="white", label=unique_functionals[i])
        r=[x + barWidth for x in r]
        i+=1
    axs[1].set_xticks([r + midFac*barWidth for r in range(len(e_coh[0]))])
    axs[1].set_xticklabels(unique_species)
    axs[1].set_ylabel("Rel. error in $E_{coh}\, (\%)$")


    #colors = iter(cm.rainbow(linspace(0, 1, len_func)))

    #r = arange(len(b_0[0]))
    #i=0
    #for b_sub in b_0:
    #    #col = next(colors)
    #    col = colors[i]
    #    axs[2].bar(r, b_sub, color=col, width=barWidth, edgecolor="white", label=unique_functionals[i])
    #    r=[x + barWidth for x in r]
    #    i+=1
    #axs[2].set_xticks([r + midFac*barWidth for r in range(len(b_0[0]))])
    #axs[2].set_xticklabels(unique_species)
    #axs[2].set_ylabel("Rel. error in $B_0$\, (\%)")


    #plt.legend(bbox_to_anchor=(0, 3.5),loc='upper left')
    plt.legend()
    #plt.show()
    fig.savefig(name,format="pdf")
    plt.close(name)
    print('Plot ',name,' was created.\n')

    return 

#def data_plot_twoD(param_list,v_list,e_list,lm_poly, *fmt_list):
#def data_plot_twoD(result_list):
def data_plot_twoD(res_dict):
    """
    A helper function to make a graph 

    Parameters
    ----------
    param_list : Axes
        The axes to draw to

    v_list : array
       The x data

    e_list : array
       The y data

    *fmt_list : dict
       Dictionary of kwargs to pass to ax.plot

    Returns
    -------
    out : list
        list of artists added
    """
    fig_width=25
    fig_height=20

    #v_list=result_list[-1];    e_list=result_list[-4];    param_list=result_list[-3]; lm_poly=result_list[2];BM_poly=result_list[1];
    v_list=res_dict['v_list']; e_list=res_dict['e_list']; param_list=res_dict['param_list']; lm_poly=res_dict['poly'];BM_poly=res_dict["BM_poly"]

    matplotlib.rcParams.update(matplotlib.rcParamsDefault)

    # FOR LATEX PLOTTING
    matplotlib.rcParams.update({'text.usetex': True,
                      'axes.labelsize': 18, # fontsize for x and y labels (was 11)
                      'axes.titlesize': 18,
                      #'text.fontsize': 8, # was 10
                      'legend.fontsize': 14, # was 10
                      'xtick.labelsize': 18,
                      'ytick.labelsize': 18,
                      'figure.figsize': [fig_width,fig_height],
                      'font.family': 'serif'})

    index_emin = e_list.index(min(e_list))

    fig = plt.figure(figsize=(fig_width, fig_height))
    
    ext_width  = 2*6.9
    ext_height = ext_width*(sqrt(5)-1)/2
    #fig2 = plt.figure(figsize=(ext_width,ext_height))
    fig3 = plt.figure(figsize=(ext_width,ext_height))
    fig.suptitle("Data and volume plot")

    # FIRST SUBPLOT
    ax      = fig.add_subplot(3,1,1)
    ax.title.set_text('Calculation results')

    ax.tricontour(param_list[0], param_list[1], e_list, levels=14, linewidths=0.5, colors='k')
    contour = ax.tricontourf(param_list[0], param_list[1], e_list, levels=14, cmap='RdBu_r')
    pts = ax.scatter(param_list[0], param_list[1],color=[0.0,0.0,0.0,1],s=10)

    fig.colorbar(contour,ax=ax)
    ax.set_xlabel(r"a [Å]")
    ax.set_ylabel(r"c/a [N/A]")
    ax.set_xlim(min(param_list[0])-0.02,max(param_list[0])+0.02)
    ax.set_ylim(min(param_list[1])-0.02,max(param_list[1])+0.02)
    nbr_unique_y = len(unique([round(e,2) for e in param_list[1]]))
    nbr_unique_x = len(unique([round(e,2) for e in param_list[0]]))
    ax.set_yticks(linspace(min(param_list[1]),max(param_list[1]),nbr_unique_y))
    ax.set_xticks(linspace(min(param_list[0]),max(param_list[0]),nbr_unique_x))
    # SELECTOR FOR FIRST PLOT
    selector = SelectFromCollection(ax, pts)

    def accept(event,e_list):
        if event.key == "enter":
            print("Selected points:")
            print(selector.xys[selector.ind])
            selector.disconnect()
            ax.set_title("")
            fig.canvas.draw()
            result_list_two=two_dim_bulkmod("replot","re_func",[array(param_list[0])[selector.ind],array(param_list[1])[selector.ind]],array(v_list)[selector.ind],array(e_list)[selector.ind],[0],"hcp")
            print() 
            print("\t\t Fourth Order Poly. Fit")
            print_ascii_output_table(result_list_two[3],"hcp","poly")
            print()
            return 

    cid = fig.canvas.mpl_connect("key_press_event", lambda event: accept(event,e_list))
    ax.set_title("Calculation results.\n Press enter to accept selected points.")

    # SECOND SUBPLOT *AS STAND-ALONE PLOT*
    #ax      = fig2.add_subplot(1,1,1)
    #ax.title.set_text('Energy vs a,c')
    #ax.tricontour(param_list[0], param_list[1], lm_poly.eval(x=[param_list[0],param_list[1]]), levels=14, linewidths=0.5, colors='k')
    #contour = ax.tricontourf(param_list[0], param_list[1], lm_poly.best_fit, levels=14, cmap='RdBu_r')
    #ax.plot(param_list[0], param_list[1],'ko',ms=3)
    #fig.colorbar(contour,ax=ax)
    #ax.set_xlabel(r"a [Å]")
    #ax.set_ylabel(r"c/a [N/A]")

    # SECOND SUBPLOT *AS STAND-ALONE PLOT*
    ax      = fig.add_subplot(3,1,2)
    ax.title.set_text('Polynomial on calc set')
    ax.tricontour(param_list[0], param_list[1], lm_poly.eval(x=[param_list[0],param_list[1]]), levels=14, linewidths=0.5, colors='k')
    contour = ax.tricontourf(param_list[0], param_list[1], lm_poly.best_fit, levels=14, cmap='RdBu_r')
    ax.plot(param_list[0], param_list[1],'ko',ms=3)
    fig.colorbar(contour,ax=ax)
    ax.set_xlabel(r"a [Å]")
    ax.set_ylabel(r"c/a [N/A]")

    # THIRD SUBPLOT
    #BM_poly = bmfit_1D(v_list,e_list)
    #v_grid  = linspace(v_list[0],v_list[-1],1000)
    #ax      = fig.add_subplot(2,1,4)
    #ax.title.set_text('Murnaghan (volume) fit')
    #ax.plot(v_list,e_list,"kx",label="Data")
    #ax.plot(v_grid,ameos(v_grid,BM_poly[0]),"-",label="Murnaghan fit")
    #ax.set_xlabel(r"V  [Å$\displaystyle ^3$]")
    #ax.set_ylabel(r"Energy [eV]")
    #ax.legend(loc="best")

    # vectors for eval of poly
    y = linspace(min(param_list[1]),max(param_list[1]))
    x = linspace(min(param_list[0]),max(param_list[0]))
    xv, yv = meshgrid(x,y)
    z=lm_poly.eval(x=[xv,yv])

    # 3D plot
    ax      = fig3.add_subplot(111, projection="3d")
    ax.title.set_text('Energy vs vol')
    ax.plot_surface(xv,yv,z, rstride=1, cstride=1, cmap='RdBu_r', edgecolor='none')
    ax.set_xlabel(r"a [Å]")
    ax.set_ylabel(r"c/a [N/A]")

    # FOURTH SUBPLOT
    ax      = fig.add_subplot(3,1,3)
    ax.title.set_text('Polynomial on evenly spaced grid')
    ax.contour(xv, yv, z, levels=14, linewidths=0.5, colors='k')
    contour  = ax.contourf(x, y, z, levels=14, cmap='RdBu_r')
    ax.plot(param_list[0], param_list[1],'ko',ms=3)

    fig.colorbar(contour,ax=ax)
    ax.set_xlabel(r"a [Å]")
    ax.set_ylabel(r"c/a [N/A]")
    return fig

def pdf_plot_twoD(plot_name,res_dict):
    """
    A helper function to make a graph 

    Parameters
    ----------
    param_list : Axes
        The axes to draw to

    v_list : array
       The x data

    e_list : array
       The y data

    *fmt_list : dict
       Dictionary of kwargs to pass to ax.plot

    Returns
    -------
    out : list
        list of artists added
    """
    fig_width=15
    fig_height=10

    #v_list=result_list[-1];    e_list=result_list[-4];    param_list=result_list[-3]; lm_poly=result_list[2];BM_poly=result_list[1];
    v_list=res_dict['v_list']; e_list=res_dict['e_list']; param_list=res_dict['param_list']; lm_poly=res_dict['poly'];BM_poly=res_dict["BM_poly"]; geo=res_dict["geo"]

    matplotlib.rcParams.update(matplotlib.rcParamsDefault)

    # FOR LATEX PLOTTING
    matplotlib.rcParams.update({'text.usetex': True,
                      'axes.labelsize': 22, # fontsize for x and y labels (was 11)
                      'axes.titlesize': 22,
                      #'text.fontsize': 8, # was 10
                      'legend.fontsize': 18, # was 10
                      'xtick.labelsize': 20,
                      'ytick.labelsize': 20,
                      'figure.figsize': [fig_width,fig_height],
                      'font.family': 'serif'})

    index_emin = e_list.index(min(e_list))

    fig = plt.figure('Polynomial evaluated on calculation grid',figsize=(fig_width, fig_height))
    fig2 = plt.figure('Calculation results',figsize=(fig_width, fig_height))
    
    # PLOT poly on grid
    ax      = fig.add_subplot(1,1,1)
    ax.tricontour(param_list[0], param_list[1], lm_poly.eval(x=[param_list[0],param_list[1]]), levels=14, linewidths=0.5, colors='k')
    contour = ax.tricontourf(param_list[0], param_list[1], lm_poly.best_fit, levels=14, cmap='RdBu_r')
    ax.plot(param_list[0], param_list[1],'ko',ms=3)
    fig.colorbar(contour,ax=ax)
    #nbr_unique_y = len(unique([round(e,2) for e in param_list[1]]))
    #nbr_unique_x = len(unique([round(e,2) for e in param_list[0]]))
    #ax.set_yticks(linspace(min(param_list[1]),max(param_list[1]),nbr_unique_y))
    #ax.set_xticks(linspace(min(param_list[0]),max(param_list[0]),nbr_unique_x))
    ax.set_xlabel(r"a [Å]")
    if (geo=="3L-ab"):
        ax.set_ylabel(r"b [Å]")
    elif (geo=="2L"):
        ax.set_ylabel(r"c/a [N/A]")

    # PLOT data on grid
    ax      = fig2.add_subplot(1,1,1)
    ax.tricontour(param_list[0], param_list[1], e_list, levels=14, linewidths=0.5, colors='k')
    contour = ax.tricontourf(param_list[0], param_list[1], e_list, levels=14, cmap='RdBu_r')
    pts = ax.scatter(param_list[0], param_list[1],color=[0.0,0.0,0.0,1],s=10)
    fig.colorbar(contour,ax=ax)
    ax.set_xlabel(r"a [Å]")
    if (geo=="3L-ab"):
        ax.set_ylabel(r"b [Å]")
    elif (geo=="2L"):
        ax.set_ylabel(r"c/a [N/A]")

    plt.show()
    name = '%s_poly.pdf'%(plot_name)
    name2 = '%s_calc.pdf'%(plot_name)
    fig.savefig(name,format="pdf")
    fig2.savefig(name2,format="pdf")
    plt.close(name)
    plt.close(name2)
    print('Plot ',name,' was created.\n')
    print('Plot ',name2,' was created.\n')

    return fig
