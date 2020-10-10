from   scipy import *
from   scipy.optimize import leastsq,curve_fit
from   lmfit import Model
from numpy import argmin

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
            1*( del1111*(x[0]-a0)**4 + del2222*(x[1]-c0)**4 + 1*( del1112*(x[0]-a0)**3*(x[1]-c0)+ del1222*(x[0]-a0)*(x[1]-c0)**3 ) +1*del1122*(x[0]-a0)**2*(x[1]-c0)**2 )


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
