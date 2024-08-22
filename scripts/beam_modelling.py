import numpy as np
import healpy as hp
from scipy.special import jv

def TaperAngle(X, a, b, c, d, f):
    return a + b*X + c*X**2 + d/X + f/X**2
    
def pseudo_BINGO_FWHM(nu, nu_ref=1100, fwhm_ref=40, out_unit="arcmin"):
    cte = np.array([-1.46802646e+02,  1.10877204e-01, -2.96922720e-05,  9.04792596e+04, -1.19053639e+07])
    inc = fwhm_ref-TaperAngle(nu_ref,*cte)    
    _fwhm_ = 2*(TaperAngle(nu, *cte) + inc)
    if out_unit=="arcmin": return _fwhm_
    elif out_unit=="deg": return _fwhm_/60
    elif out_unit=="radians": return np.radians(_fwhm_/60)
        
def jinc(x):
    j = np.zeros_like(x)
    ind0    = np.where(x==0)[0]
    indnot0 = np.where(x!=0)[0]    
    if ind0.size:
        j[ind0]=0.5
    j[indnot0] = jv(1, x[indnot0])/x[indnot0]
    return j
    
#[Warning!] Matshawule calculated his polynomial coefficients in MHz and sine parameters with Hz
def fwhm_cosine_coef_model(nu=None, type_=None):# Matshawule et al (2021), Eq. (2)
    an  = np.array([ 3.40234907e-21, -3.02516490e-17,  1.17019761e-13, -2.57168633e-10,
                     3.51130410e-07, -3.04953935e-04,  1.64488782e-01, -5.03702034e+01,6.70428133e+03])
    coefsum = np.poly1d(an)(nu)
    A   = 0.1  #[arc min]
    T   = 20   #[MHz]
    if type_=="ripple":
        return coefsum+(A/60)*np.sin(2*np.pi*nu/T)
    if  type_=="smooth":
        return coefsum
    if  type_=="constant":
        return 1.16*np.ones_like(nu)
    
def fwhm_modelling(nu=None,type_=None, D=13.5, in_degree=True, if_fwhm_fixed=None):
    c    = 299792458 #m s-1
    #nu  = nu*1e6
    lbda = c/(nu*1e6)
    if in_degree:
        lbda_D = np.degrees(lbda/D)
    else:
        lbda_D = lbda/D
    if type_=="irfan2023":
        fwhm = 1.2*(1280/nu) #in degree
        if not in_degree:
            fwhm = np.radians(fwhm)        
    elif type_=="ripple" or type_=="smooth" or type_=="constant":
        fwhm = fwhm_cosine_coef_model(nu, type_)*lbda_D
    elif type_=="carucci2020":
        fwhm = lbda_D
    elif type_=='fixed':
        fwhm = if_fwhm_fixed
    else:
        fwhm = 1.22*lbda_D
    return fwhm
    
def beam_function(type_="gaussian", fwhm=None, theta_=None): 
    if   type_=="gaussian":
        b  = (theta_/fwhm)**2
        return np.exp(-4*np.log(2)*b)
    elif type_=="jinc":
        theta_ = theta_/fwhm
        return 4*jinc(4*theta_)**2
    elif type_=="cosine":
        theta_ = 1.189*theta_
        theta_ = theta_/fwhm
        c = np.cos(np.pi*theta_)
        c = c/(1-4*theta_**2)
        return c**2

def bl_function(type_="gaussian", fwhm=None, lmax=None, theta_=None, nu=1100, input_unit="radian"): 
    if input_unit=="radian":
        pass
    elif input_unit=="degree":
        fwhm   = np.radians(fwhm)
        theta_ = np.radians(theta_)
    elif input_unit=="arcmin":
        fwhm = np.radians(fwhm/60)
        theta_ = np.radians(theta_/60)
    else:
        pass   
    bl  = beam_function(type_,fwhm=fwhm, theta_=theta_)
    bl  = hp.beam2bl(bl,theta_,lmax)    
    return bl/bl.max()
