"""
ml_read_data module
written by Lucia Klarmann
contains functions that read in
   outputfiles from MCMax,
   modified .sav from Koen Maaskant
returns read data in a (hopefully) usefull format
"""

import numpy as np
#import matplotlib as mpl
import scipy
#from scipy import constants
from astropy.io import fits
import pandas as pd
from scipy import constants as const


def GHz2my(f):
    """
    convert frequncy into wavelength

    input: Frequency in GHZ, float or np array
    output: wavelength in micron, float or np array
    """
    f_SI = f*1e9
    w_SI = const.c/f_SI
    w = w_SI * 1e6

    return (w)


def Jy2W(Fd, f):
    """
    convert flux density to flux
    input: Flux density in Jy
           wavelength in micron
    output: Flux in W/m**2
    """
    f_SI = f*1e9
    F = Fd * 10**(-26) * f_SI

    return(F)


def Jy2erg(Fd, f):
    """
    convert flux density to flux
    input: Flux density in Jy
           wavelength in micron
    output: Flux in erg/m**2/s/my    
    """

    F = Jy2W(Fd, f)
    F = F * 1e3  # my, not meter!!!
    return (F)


def raytrace_data(path, name):
    """
    reads and returns raytraced data from MCMax

    keyword arguments:
    path -- path to spectrum??.?.dat file
    name -- spectrum??.?.dat

    returns:
    wavelenght lambda [my]
    flux density lambda [erg s^-1 cm^-2 my^-1]
    flux lambda [erg s^-1 cm^-2]
    flux from scatter lambda [erg s^-1 cm^-2]
    """
    data = np.loadtxt(path + name)  # read in spectrum

    wavelength = data[:,0]

    flux_density_nu = data[:,1]
    flux_density_lambda = (flux_density_nu *
                           1e-23 * constants.c * 1e-2 /
                           (wavelength * wavelength * 1e-8))

    flux_density_fromscatter_nu = data[:,2]
    flux_density_fromscatter_lambda = (flux_density_fromscatter_nu *
                                       1e-25 * const.c /
                                       (wavelength * wavelength * 1e-8))

    flux_lambda = flux_density_lambda * wavelength
    flux_fromscatter_lambda = flux_density_fromscatter_lambda * wavelength

    return (wavelength, flux_density_lambda, flux_lambda,
            flux_fromscatter_lambda)


def MCSpec_data(path, angle):
    """
    reads and returns MCSpec data from MCMax

    keyword arguments:
    path -- path to MCSpec??.?.dat file
    select_list -- numbers of files that are to be plotted

    returns:
    wavelenght lambda [my]
    flux density lambda [erg s^-1 cm^-2 my^-1]
    flux lambda [erg s^-1 cm^-2]
    #flux from scatter lambda [erg s^-1 cm^-2]
    """
    data = np.loadtxt(path + 'MCSpec' + str(angle) + '.dat')  # read in spectrum

    wavelength = data[:,0]
    flux_density_nu = data[:,1]
    flux_density_lambda = (flux_density_nu * 1e-23 *
                           const.c * 1e-2 / (wavelength * wavelength * 1e-8))
#    flux_density_fromscatter_nu = data[:,2]
#    flux_density_fromscatter_lambda = (flux_density_fromscatter_nu * 1e-23 *
#                                       constants.c / (wavelength * wavelength * 1e-8))
    flux_lambda = flux_density_lambda * wavelength
#    flux_fromscatter_lambda = flux_density_fromscatter_lambda * wavelength

    return (wavelength, flux_density_lambda, flux_lambda)#, flux_fromscatter_lambda)

def MCSpec_star_scatter(path, angle):
    """
    reads and returns MCSpec data from MCMax

    keyword arguments:
    path -- path to MCSpec??.?.dat file
    select_list -- numbers of files that are to be plotted

    returns:
    wavelenght lambda [my]
    flux density lambda [erg s^-1 cm^-2 my^-1]
    flux lambda [erg s^-1 cm^-2]
    #flux from scatter lambda [erg s^-1 cm^-2]
    """
    data = np.loadtxt(path + 'MCSpec' + str(angle) + '.dat')  # read in spectrum

    wavelength = data[:,0]
    flux_density_star = data[:,5]
    flux_density_scatter = data[:,3]
    flux_density_lambda_star = (flux_density_star * 1e-23 *
                           const.c * 1e-2 / (wavelength * wavelength * 1e-8))
    flux_density_lambda_scatter = (flux_density_scatter * 1e-23 *
                           const.c * 1e-2 / (wavelength * wavelength * 1e-8))
#    flux_density_fromscatter_nu = data[:,2]
#    flux_density_fromscatter_lambda = (flux_density_fromscatter_nu * 1e-23 *
#                                       constants.c / (wavelength * wavelength * 1e-8))
    flux_lambda_star = flux_density_lambda_star * wavelength
    flux_lambda_scatter = flux_density_lambda_scatter * wavelength    
#    flux_fromscatter_lambda = flux_density_fromscatter_lambda * wavelength

    return (wavelength, flux_lambda_star, flux_lambda_scatter)#, flux_fromscatter_lambda)

# def MCSpec_all_data_NIR(path, angle):
#     """
#     Reads the MCSpec of a certain angle and gives back Fluxes in the NIR

#     keyword arguments:
#     path -- path to MCSpec??.?.dat file
#     angle -- inclination angle of MCSpec??.?.dat (the ?)

#     returns:
#     wavelength wave [my]
#     total flux flux_lambda
#     stellar flux flux_star
#     scattered flux flux_scat
#     """
#     data = np.loadtxt(path + 'MCSpec' + str(angle) + '.dat')  # read in spectrum

#     wavelength = data[:,0]
#     flux_density_star = data[:,5]
#     flux_density_scatter = data[:,3]
#     flux_density_lambda_star = (flux_density_star * 1e-23 *
#                            constants.c * 1e-2 / (wavelength * wavelength * 1e-8))
#     flux_density_lambda_scatter = (flux_density_scatter * 1e-23 *
#                            constants.c * 1e-2 / (wavelength * wavelength * 1e-8))

#     flux_lambda_star = flux_density_lambda_star * wavelength
#     flux_lambda_scatter = flux_density_lambda_scatter * wavelength    
#     flux_fromscatter_lambda = flux_density_fromscatter_lambda * wavelength

#     flux_rel_star_NIR = np.asarray(flux_lambda_star[58])/np.asarray(fl)

#     return()

def obs_data(path, name, method, list_methods=None):
    """
    Reads and returns the modified data from Koen.

    keyword arguments:
    path -- path to .sav file
    name -- name of .sav file
    method -- which data to give back:
    method_list -- if true: prints list of observed data options

    Please note that kurucz data has no errors and returns a 2-Tuple!

    returns:
    wavelength [my]
    flux [erg s^-1 cm^-2]
    flux error [erg s^-1 cm^-2]
    for chosen observation method
    """
    method_list = ['kurucz', 'kurucz_unred',
                   'phot', 'phot_other', 'phot_akari', 'phot_pacs',
                   'spitzer', 'spitzer_rebinned']
    raw_obs_dict = scipy.io.readsav(path + name, python_dict=True)

    ## #print keys if method list is true
    ## if method_list:
    ##     for key, value in raw_obs_dict.iteritems():
    ##         print(key)

    # print available methods if method list
    if list_methods:
        print(method_list)

    if method not in method_list:
        print('please use any word as fourth argument to get list of methods,'
              'this method is not defined')
        return

    # no errors for kurucz data
    pattern = 'kurucz'

    if pattern in method:
        wave = raw_obs_dict[method+'_wave']
        flux = raw_obs_dict[method+'_flux']
        return (wave, flux)

    # for all other observation methods:
    wave = raw_obs_dict[method+'_wave']
    flux = raw_obs_dict[method+'_flux']
    err = raw_obs_dict[method+'_err']
    return (wave, flux, err)


def scaleheights_Q_data(path, name, unit='r'):
    """
    Reads scalheight.dat from MCMax and returns data

    keyword arguments:
    path -- path to scaleheight.dat file
    name -- name of scaleheight.dat  file
    unit -- unit of scaleheigt, standard is r, set unit='AU' for AU

    returns:
    radius [AU]
    scaleheight gas [r or AU]
    scaleheight dust [r or AU]
    Toomre Q
    """
    data = np.loadtxt(path + name)  # read in data

    radius = data[:,0]
    h_gas = data[:,1]
    h_dust = data[:,2]
    Q_t = data[:,3]

    if unit == 'AU':
        return (radius, h_gas, h_dust, Q_t)

    elif unit == 'r':
        h_gas = h_gas / radius
        h_dust = h_dust / radius
        return (radius, h_gas, h_dust, Q_t)

    else:
        print('The units you choose can not be returned.')
        print('For scaleheights in AU, use unit="AU"')
        print('For scaleheights normed to radius use unit="r"'
              'or omit the unit keyword')
        return


def selected_kappa_data(path, name, select_list=[],
                        label_path='', label_name=''):
    """
    reads kappas.dat and returns selected kappas

    keyword arguments:
    path -- path to kappas.dat file
    name -- name of kappas.dat file
    select_list -- numbers of kappas to be read, ommit to read all

    returns
    wave --  wavelength [my]
    kappas_scat -- scatter opacity [cm^2 g^-1]
    kappas_abs -- absorbtion opacity [cm^2 g^-1]
    label -- grain sizes for each lambda for legend u.ae.

    TODO: Make mrn file optional, it is not always created
    Do not use selct list. Easier to get all and select while plotting!
    """

    if not label_path:
        label_path = path
        #print ('da')

    if not label_name:
        label_name = 'mrn.dat'

    print 'here'
    print path+name
    data = np.genfromtxt(path + name, float)
    label = np.genfromtxt(label_path + label_name)

    wave = data[:,0]
    kappas2 = len(data[0])-1
    number_kappas = kappas2 / 2

    if not select_list:
        
        #print('2=kappas2 %i' % kappas2)  # for debugging
        #print('number_kappas %i' % number_kappas)  # for debugging
        kappas_abs = data[:, range(1, number_kappas+1)]
        kappas_scat = data[:, range(number_kappas+1, kappas2+1)]
        #print (label)

        return (wave, kappas_abs, kappas_scat, label)

    else:
        kappas_abs = data[:,select_list]
        select_list = np.array(select_list)
        select_list = select_list + number_kappas
        kappas_scat = data[:,select_list]
        #print (kappas_abs)
        #print (kappas_scat)
        #return
        return(wave, kappas_abs, kappas_scat, label)


#####################################################################

def kappa_data_nolabel(path, name):
    """
    reads kappas.dat and returns selected kappas

    keyword arguments:
    path -- path to kappas.dat file
    name -- name of kappas.dat file
    select_list -- numbers of kappas to be read, ommit to read all

    returns
    wave --  wavelength [my]
    kappas_scat -- scatter opacity [cm^2 g^-1]
    kappas_abs -- absorbtion opacity [cm^2 g^-1]
    label -- grain sizes for each lambda for legend u.ae.

    TODO: Make mrn file optional, it is not always created
    Do not use selct list. Easier to get all and select while plotting!
    """
    print path+name
    data = np.genfromtxt(path + name, float)


    wave = data[:,0]
    kappas2 = len(data[0])-1
    number_kappas = kappas2 / 2


        
    #print('2=kappas2 %i' % kappas2)  # for debugging
    #print('number_kappas %i' % number_kappas)  # for debugging
    kappas_abs = data[:, range(1, number_kappas+1)]
    kappas_scat = data[:, range(number_kappas+1, kappas2+1)]
    #print (label)

    return (wave, kappas_abs, kappas_scat)



####################################################################


def denstemp_data(path, name, unit='r', mesh='mesh'):
    """
    read density and temperature from denstemp.fits and return it

    keyword arguments:
    path -- path to denstemp.fits file
    name -- name of denstemp.fits file
    unit -- unit of angle/z coordinate, standard is z [r]

    returns
    radius [AU]
    height coordinate -- angle [rad], z [AU] or [r]
    density [g/cm^-3]
    temperature [K]

    This data spans a large range of values. Use log for plots.
    DONE integrate mesh option for fast plot
    Please note that the coordinates are returned as a meshgrid
    as long as the mesh keyword is not set to 'no'
    """

    data_list = fits.open(path + name)

    radius = data_list[0].data[0][0,:]  # read in radius
    angle = data_list[0].data[1][:,1]  # declination starts at pi/2!

    dens = data_list[1].data  # read in density
    temp = data_list[2].data  # read in temperatur
    print(radius.shape)
    print(angle.shape)

    if mesh == 'mesh':
        radius, angle = np.meshgrid(radius, angle)
        print(radius.shape)
        print(angle.shape)

    elif mesh == 'no':
        pass

    else:
        print('The coordinates can not be returned in that format')
        print('please use mesh="no" for 1d arrays')
        print('omit or use mesh="mesh" for np.meshgrid')

    if unit == 'rad':
        return(radius, angle, dens, temp)

    elif unit == 'r':
        z_over_r = 1. / np.tan(angle)
        return(radius, z_over_r, dens, temp)

    elif unit == 'AU':
        z = (1. / np.tan(angle)) * radius
        return(radius, z, dens, temp)

    else:
        print('the units you have choosen can not be returned')
        print('please use unit = "rad" for declination in rad')
        print('use unit = "AU" for z in AU and unit = "r" or omit'
              'keyword for z in r')
        return


def height_data(path, name, select='TODO', set_zero='no', unit='r'):
    """
    read perpendicular optical depth from height.dat and return

    keyword arguments:
    path -- path to height.dat file
    name -- name of height.dat file
    select -- creat tau_X_X with MCMax and select it here. TODO
    set_zero -- set negative values zero for plotting
    unit -- of z coordinate where tau = 1, choose "AU" or standard [r]

    returns
    radius [AU]
    height for tau= 1 in [AU] or [r] for wavelength 0.55 or
    DONE selected lambda: use multi_height_data
    DONE ask Michiel about tauR in last column
    """
    data = np.loadtxt(path + name)  # read in data

    radius = data[:,0]
    tau_0_55 = data[:,1]

    if unit == 'AU':
        pass

    elif unit == 'r':
        tau_0_55 = tau_0_55 / radius

    else:
        print('the units you have choosen can not be returned')
        print('please use unit = "AU" for z in AU and unit = "r" or omit'
              'keyword for z in r')
        return
    
    if set_zero == 'no':
        return (radius, tau_0_55)

    elif set_zero == 'set_zero':
        tau_0_55[np.where(tau_0_55 < 0)] = 1e-34
        return (radius, tau_0_55)

    else:
        print('Please use set_zero = "set_zero" to set negative'
              'values of tau zero')
        print('Omit keyword or use set_zero = "no" for original data')

def multi_height_data(path, name, number, set_zero='no', unit='r'):
    """
    read perpendicular optical depth from height.dat and return

    keyword arguments:
    path -- path to height.dat file
    name -- name of height.dat file
    number - number of wavelenthes for which MCMax created tau=1
    set_zero - set values <0 zero for plotting or log
    unit -- of z coordinate where tau = 1, choose "AU" or standard [r]

    returns
    radius [AU]
    height for tau= 1 in [AU] or [r] for wavelength 0.55 or
    """
    data = np.loadtxt(path + name)

    radius = data[:,0]
    tau1 = []
    #tau1 = np.asarray(tau1)
    for i in xrange(number):
        j=i+1
        tau1.append(data[:,j])

        if unit == 'AU':
            pass
        elif unit == 'r':
 #          # carful: i indices lists, not elemnts
            # i=0 includes the intere first tau column of height.dat
            tau1[i] = tau1[i]/radius
        else:
            print('the units you have choosen can not be returned')
            print('please use unit = "AU" for z in AU and unit = "r" or omit'
              'keyword for z in r')
            return
        
        if set_zero == 'set_zero':
            tau1[i,:][np.where(tau1[i,:] < 0)] = 0
            
        
    return(radius, tau1)

def heightR_data(path, name, select='TODO', unit='r'):
    """
    read radial optical depth from heightR.dat and return

    keyword arguments:
    path -- path to heightR.dat file
    name -- name of heightR.dat file
    select -- creat tau_X_X with MCMax and select it here. TODO
    unit -- of z coordinate where tau = 1, choose "AU" or standard [r]

    returns
    x-coordinate (r) [AU]
    height for tau= 1 in [AU] or [r] for wavelength 0.55 or
    TODO selected lambda
    """
    data = np.loadtxt(path + name)  # read in data

    tau_x = data[:,0]
    tau_z = data[:,1]

    if unit == 'AU':
        pass

    elif unit == 'r':
        tau_z = tau_z / tau_x

    else:
        print('the units you have choosen can not be returned')
        print('please use unit = "AU" for z in AU and unit = "r" or omit'
              'keyword for z in r')
        return

    return (tau_x, tau_z)




def multi_heightR_data(path, name, number, unit='r'):
    """
    return radial optical depth of selcted wavelenth from heightR.dat

    keyword arguments:
    path - path to heightR.dat file
    name - name of heightR.dat file
    number - number of wavelengthes, defined in .in file MCMax 
    unit - of z coordinate where tau=1, choose "AU" or standard [r]

    returns
    x-coordinate in [AU] or [r], same for each wavelength
    height for tau=1 in [AU] or [r], same for each wavelength
    """
    data = np.loadtxt(path + name)  # read in data

    tau_x1 = data[:,0]
    tau_z1 = data[:,1]

    if unit == 'AU':
        pass

    elif unit == 'r':
        tau_z1 = tau_z1 / tau_x1

    else:
        print('the units you have choosen can not be returned')
        print('please use unit = "AU" for z in AU and unit = "r" or omit'
          'keyword for z in r')
        return
    
    if number == 1:
        return (tau_x1, tau_z1)
        
    else:
        tau_x = []
        tau_z = []
        for i in xrange(number-1):
            j = 2*i+3
            tau_x.append(data[:,j])
            tau_z.append(data[:,j+1])
            

            if unit == 'AU':
                pass
            elif unit == 'r':
                tau_z[i]= tau_z[i] / tau_x[i]

            else:
                print('the units you have choosen can not be returned')
                print('please use unit = "AU" for z in AU and unit = "r" or omit'
                      'keyword for z in r')
                return

    return(tau_x1, tau_z1, tau_x, tau_z)
    
def surfacedens_data(path, name):
    """
    reads surfacedensity from surfacedens.dat and returns it

    keyword arguments:
    path -- path to heightR.dat file
    name -- name of heightR.dat file

    radius [AU]
    surface density [g/cm^2]

    """
    data = np.loadtxt(path + name)  # read in data

    radius = data[:,0]
    surfacedens = data[:,1]

    return (radius, surfacedens)


def radtau_data(path, name):
    """
    reads stellar and 0.55 optical depth of mp from radtau.dat

    keyword arguments:
    path -- path to radtau.dat file
    name -- name of radtau.dat file

    radius [AU]
    tau [unitless]
    """
    data = np.loadtxt(path + name)

    radius = data[:,0]
    tau_mp_stellar = data[:,1]
    tau_mp_0_55mu = data[:,2]

    return (radius, tau_mp_stellar, tau_mp_0_55mu)

def radtau_data_multi(path, name, number):
    """
    reads the tau in midplane for the lam0n wavelength from the .in file

    keyword argumetns:
    path -- path to radtau.dat file
    name -- nam eof radtau.dat file

    radius [AU]
    tau [unitless] for each wavelength lam0n
    """
    #print path
    #print name
    #print number
    
    data = np.loadtxt(path+name)
    #print data

    radius = data[:,0]
    
    tau = []
    for i in xrange(number):
        tau.append(data[:,i+1])

    return(radius, tau)
        

#input files for MCMax from Paola, from Tils Code

def radgrid_data(path, name):
    """
    reads in the radial grid points from radgrid.dat

    keyword arguments:
    path -- path to radgrid_XXX.dat file
    name -- name of radgrid_XXX.dat file

    radius [AU]
    """

    data = np.loadtxt(path + name)
    radius = data[:]


    return (radius)
        

def raddens_paola(path, name):
    """
    reads in Paola's raddens files

    keyword arguments:
    path -- path to raddens_XXX.dat file
    name -- name of raddens_XXX.dat file

    radius [AU]
    surface density [g/cm^2]
    """
    data = np.loadtxt(path + name)

    radius = data[:,0]
    surfacedens_data = data[:,1]

    return (radius, surfacedens_data)


def surfdens_paola(path, name):
    """
    reads in Paola's surfdens files, surfdenses for each particle

    keyword arguments:
    path -- path to surfdens_XXX.dat file
    name -- name of surfdens_XXX.dat file

    radius [AU]
    surface density [g/cm^2]
    """

    data = np.loadtxt(path + name)

    #number_dust_part = len(data[0,:])-2
    dust_parts_dens = data[:,1:-1]

    radius = data[:,0]
    gasdens = data[:,-1]

    return (radius, gasdens, dust_parts_dens)
        
def basevis_pionier(path, name, number_wave):
    """
    reads in the interferometry pionier data, baseline and visibilities
    always use three wavelengths

    keyword arguments:
    path -- path to phasevisXX.XPIONIER.dat file
    name -- name of phasevisXX.XPIONIER.dat file
    number_wave -- number of wavlengthes traced for file

    baseline []
    full disk at last wavelength
    visibilities for all wavelength
    phase for all wavelength
    """

    data = np.loadtxt(path + name)
    #2:4
    base = data[:,0]
    full_disk = data[:,1]
    vis = data[:, 2:2+number_wave]
    phase = data[:, 2+number_wave:2+2*number_wave]

    return(base, full_disk, vis, phase)


def vis_matisse(path, name, number_basexangle):
    """
    reads in the interferometry pionier data, baseline and visibilities
    always use three wavelengths

    keyword arguments:
    path -- path to visXX.XMATSPEC_X.dat file
    name -- name of visXX.XMATSPEC_X.dat file
    number_basexangle -- number of baselines times position angles

    wavelength []
    full disk at last wavelength
    visibilities for all baselines, angels
    phase for all baslines, angles
    """

    data = np.loadtxt(path + name)
    #2:4
    base = data[:,0]
    full_disk = data[:,1]
    vis = data[:, 2:2+number_basexangle]
    phase = data[:, 2+number_basexangle:2+2*number_basexangle]

    return(base, full_disk, vis, phase)

def spitzer_noerr(path, name):
    """
    gets the spitzer data extracted from idl, needs to be
    in two columns, until now only hd100453

    keyword arguments:
    path -- path to two column data file
    name -- name of datafile

    wavelength -- wavelength in micron
    flux -- flux  [erg s^-1 cm^-2] (units uncertain!) ???
    TODO: Check flux units
    """

    data = np.loadtxt(path + name)

    wave = data[:,0]
    flux = data[:,1]

    return(wave, flux)
#TODO add errors to spitzer spectrum
def spitzer(path, name):
    """
    gets the spitzer data extracted from idl, needs to be
    in two columns, until now only hd100453

    keyword arguments:
    path -- path to two column data file
    name -- name of datafile

    wavelength -- wavelength in micron
    flux -- flux  [erg s^-1 cm^-2] (units uncertain!) ???
    TODO: Check flux units
    """

    data = np.loadtxt(path + name)

    wave = data[:,0]
    flux = data[:,1]
    error = data[:,2]

    return(wave, flux, error)

def radial_flux(path, name, unit='Jy', wavelength=1.67):
    """
    reads in the radial flux distribution frrom the RImage file
    scaletype1, Jy and AU
    keyword arguments:
    
    path -- path to two column data file
    name -- name of datafile
    unit -- default in Jy
    wavelength -- default 1.67 micron

    radius -- radius [AU]
    flux -- cummulative flux  [Jy] (units uncertain!) ???
    """

    data = np.loadtxt(path + name)

    radius = data[:,0]
    flux_density = data[:,1]

    if unit == 'Jy':
        flux = flux_density
        return(radius, flux)

    elif unit == 'erg':
        flux = (flux_density * 1e-23 *
                           const.c * 1e-2 / (wavelength * wavelength * 1e-8)) *wavelength
        return(radius, flux)

def radial_flux_arcsec(path, name, unit='mJy', wavelength=1.67):
    """
    reads in the radial flux distribution frrom the RImage file
    scaletyp2, mJy and arcsec
    keyword arguments:
    
    path -- path to two column data file
    name -- name of datafile
    unit -- default in mJy
    wavelength -- default 1.67 micron

    radius -- radius [mas]
    flux -- cummulative flux  [Jy] (units uncertain!) ???
    """

    data = np.loadtxt(path + name)

    radius = data[:,0]
    flux_density = data[:,1]

    if unit == 'mJy':
        flux = flux_density
        return(radius, flux)

    elif unit == 'erg':
        flux = (flux_density * 1e-3 * 1e-23 *
                           const.c * 1e-2 / (wavelength * wavelength * 1e-8)) *wavelength
        return(radius, flux)

def read_csv_dgtau(path, name, unit='erg'):
    """
    reads in DG Tau photometry from VIZIR csv file

    arguments:
    path -- path to data file
    name -- name of data file
    unit -- unit of flux/flux density
            erg - erg/cm**2/s/my
            Watt - W/m**2
            Jy  - Jy
    returns:
    wave -- wavelength in micron
    Fd/flux flux/flux density in unit
    Fd_err/flux_err error on flux density/flux in unit
    """
    data = pd.read_csv(path+name)
    Fd = np.asarray(data['_sed_flux'])
    Fd_err = np.asarray(data['_sed_eflux'])
    f = np.asarray(data['_sed_freq'])

    wave = GHz2my(f)

    if unit=="Jy":
        return (wave, Fd, Fd_err)

    if unit=="W":
        flux = Jy2W(Fd, f)
        flux_err = Jy2W(Fd_err, f)
        return (wave, flux, flux_err)

    if unit=='erg':
        flux = Jy2erg(Fd, f)
        flux_err =  Jy2erg(Fd_err, f)
        return (wave, flux, flux_err)

    else:
        print 'please choose a correct unit. Default is erg/cm**2/s/my'
