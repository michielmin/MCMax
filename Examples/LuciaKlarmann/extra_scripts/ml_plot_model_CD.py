"""
ml_plot_model module
written by Lucia Klarmann
plots all output as single files and multifile
   outputfiles from MCMax
"""

import numpy as np
import scipy
import matplotlib as mpl
import ml_read_data as mlr
import intfits_quickread as ifq

#mpl.use('ps')
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.colors import LogNorm
#

mpl.rc('image', cmap='jet')
mpl.style.use('classic') # keeps the rainbow colors for density and temeratures


def model_multiplot(common_path, model_name, spectrum_name='spectrum35.0.dat',
                    timestep='', basevis='', PAH_free=''):
    """
    reads MCMax Output and returns one plot with eight suplots
    
    keyword arguments:
    common_path -- path to model folder
    model_name -- name of model
    spectrum_name -- name of spectrum, default is 35.0 degrees
                     set to 'no_spectrum' to get spectrum from MCSpec
    timestep -- timestep to plot from default is empty, eg. last
                computation. Takes leading zeros for one digit numbers
                for mor then 10 timesteps etc.
    basevis -- set to basevis to plot visibilities instead of inclined spectra

    returns:
    One .pdf multiplot containing
    spectrum, flux and from scatter, from star (raytrace or MCSpec)
    kappas
    temperature, inclaudes tau and tauR
    density, includes tau and tauR
    surfacedensity
    scalheights and toomre Q
    spectrum for different angles or visibilities
    zoom into tempdes, linear

    Please make sure that all files exist before calling the function
    """

    path = common_path + model_name + '/Output/'
    #Add name of star only model if needed
    model_star = 'Paola_star'
#    print 'here3'
    #plotpath_single = path
    plotpath_multi = path

    # Open plot for multiplot
    multi_fig = plt.figure(figsize=(7.5,11.0))

    # plot spectrum
    ax1 = multi_fig.add_subplot(4,2,1)
#    spectrum_name = 'spectrum35.0.dat'

    #if there has been no raytracing, use mcspec45
    if(spectrum_name == 'no_spectrum'):

        

            #plot spectrum for different inclinations

        all_angles = [10.5, 18.2, 23.6, 28.0, 31.8, 35.2, 38.4, 41.4, 44.2,
                      46.9, 49.5, 51.9, 54.3, 56.6, 58.9, 61.1, 63.3, 65.4, 67.5,
                      69.5, 71.5, 73.5, 75.5, 77.5, 79.4, 81.4, 83.3, 85.2, 87.1, 89.0]
        angle_index=[8]
        angles = [all_angles[int(i)]for i in angle_index]
    #    print(angles)

        for angle in angles:
            MCSpec_wave, MCSpec_flux_dens, MCSpec_flux = mlr.MCSpec_data(path, angle)
            ax1.plot(MCSpec_wave, MCSpec_flux, label='$\lambda F_{\lambda}$ '+str(angle)+' deg')

            MCSpec_wave2, MCSpec_flux_star, MCSpec_flux_scatter = mlr.MCSpec_star_scatter(path, angle)
            ax1.plot(MCSpec_wave2, MCSpec_flux_scatter, label='$\lambda F_{\lambda,\, \mathrm{sc}}$')
            ax1.plot(MCSpec_wave2, MCSpec_flux_star, label = '$\lambda F_{\lambda,\, \mathrm{star}}$')

        # For each long-term project, I hard code some data, e.g. photometry or spectra. Left as example, can be deleted    
            
        spitzer = False
        if spitzer == True:
            wave_spitzer, flux_spitzer = mlr.spitzer_noerr(
                '/home/lucia3/Documents/interferometry/', 'hd100453_spitzer.dat')
#            wave_spitzer, flux_spitzer = mlr.spitzer_noerr("\home\lucia3\Documents\Interferometry\", "hd100453_spitzer.dat")
            ax1.plot(wave_spitzer, flux_spitzer, '-', lw=0.5)

        photometry = False
        if photometry == True:
            wave_photo = [-0.444, -0.357, -0.260, 0.089, 0.218, 0.346, 0.576, 0.679, 1.072, 1.387]
            flux_photo = [-10.6020, -10.3750, -10.5270, -11.3150, -11.5481, -11.7045, -11.9968, -12.2599, -12.8078, -12.7717]
            wave_photo_unlog=[]
            flux_photo_unlog=[]
            wave_photo_plot=[]
            flux_photo_plot=[]
            print len(wave_photo)
            print len(flux_photo)
            for l in xrange(len(wave_photo)):
                wave_photo_unlog.append(10**wave_photo[l])
                flux_photo_unlog.append(10**flux_photo[l])
                wave_photo_plot.append(wave_photo_unlog[l])
                flux_photo_plot.append(flux_photo_unlog[l] * 1e3 * wave_photo_unlog[l])
            print wave_photo_unlog
            print flux_photo_unlog
            print wave_photo_plot
            print flux_photo_plot
            print 'still at ax1'
            ax1.plot(wave_photo_plot, flux_photo_plot, '+')

        DGTau = False
        if DGTau == True:
            wave, flux, flux_err = mlr.read_csv_dgtau('/home/klarmann/GRAVITY_GTO/DGTau/','phot.csv')
            print wave, flux
            ax1.plot(wave, flux, '+')

        HR5999 = False
        if HR5999 ==True:
            wave, flux, flux_err = mlr.read_csv_dgtau('/home/klarmann/GRAVITY_GTO/HR5999/RT/','photometry_HR5999')
            ax1.plot(wave, flux, 'x')
            
            
        l = plt.legend(prop={'size':5}, title='Spectrum', borderpad=0.5, ncol=1, loc=1)
        plt.setp(l.get_title(), fontsize='5')
        l.get_frame().set_alpha(0.5)
        plt.xscale('log')
        plt.yscale('log')
        plt.xlabel('$\lambda$ [$ \mathregular{\mu}$m]')
        plt.ylabel('$\lambda F_\lambda$ [erg cm$^{-2}$s$^{-1}$]')
        #plt.xlim([9e-2,2e3]) #DIANApaper
        #plt.ylim([2e-11,5e-8]) #DIANApaper
        #plt.xlim([1e-1,1e4])  #interferomtry      
        #plt.ylim([1e-12,1e-7]) #interferometry
        plt.xlim([1e-1, 1e4])
        plt.ylim([1e-15, 1e-7])
        plt.title(model_name)

        # plots and additional spectrum
        # here, e.g., there is and extra spectrum were the PAH emission is not included
        if PAH_free:
            (wave, flux_density_lambda,
         flux_lambda, flux_fromscatter_lambda) = mlr.raytrace_data(path, PAH_free)
            ax1.plot(wave, flux_lambda,':',  color='blue')#, wave_star, flux_lambda_star)
            ax1.plot(wave, flux_fromscatter_lambda,':',  color='green')#, wave_star, flux_lambda_star)            
            plt.plot([2e-1],[2e-3], ':', c='k', label='without PAHs')
            l = plt.legend(prop={'size':5}, title='Spectrum', borderpad=0.5, ncol=1, loc=1)
            plt.setp(l.get_title(), fontsize='5')

    # if there is a raytraiced spctrum file, this will be plotted
    else:
        (wave, flux_density_lambda,
         flux_lambda, flux_fromscatter_lambda) = mlr.raytrace_data(path, spectrum_name)

        # get star data if star is supposed to be plotted
#        star_path = common_path + model_star + "/Output/"
#        star_spectrum_name = 'spectrum15.0_star.dat'
#        (wave_star, flux_density_lambda_star,
#         flux_lambda_star, flux_fromscatter_lambda_star) = mlr.raytrace_data(star_path, star_spectrum_name)

        #fig_spec = plt.figure()
        ax1.plot(wave, flux_lambda, wave, flux_fromscatter_lambda)#, wave_star, flux_lambda_star)
    #    ax1.plot([10],[10],[10],[10])

        plt.legend(['$\lambda F_\lambda$', '$\lambda F_{\lambda,\, \mathrm{sc}}$'],
                   loc=1, borderpad=0.15, fontsize='8') #'$\lambda F_{\lambda,\, \mathrm{star}}$'],
        plt.xscale('log')
        plt.yscale('log')
        #plt.xlim([1,70])
        #plt.ylim([1e-13,5e-9])
        #plt.ylim([1e-13,5e-8])
        plt.ylim([1e-12,5e-7])

        plt.xlabel('$\lambda$ [$ \mathregular{\mu}$m]')
        plt.ylabel('$\lambda F_\lambda$ [erg cm$^{-2}$s$^{-1}$]')
        plt.title(model_name)

      # comment in this reagion to plot the kappas with automatic labelling from the mrn file
      # only useful if there are not too many grains, and they are not too large
      # carful, the mrn file uses Angstrom, not micron
#     # plot kappas with label
#     kappas_name = 'kappas.dat'
#     select_list = [0, 50, 100, 150, 170]#, 4, 8, 12, 16, 19]
#     ax2 = multi_fig.add_subplot(4, 2, 2)

#     wave, kappas_abs, kappas_scat, labeling = mlr.selected_kappa_data(path, kappas_name)

#     kappas_abs = np.transpose(kappas_abs)
#     kappas_scat = np.transpose(kappas_scat)

# #    print (wave.shape)
# #    print (kappas_abs[0].shape)
#     for i in select_list:
#         ax2.plot(wave, kappas_abs[i],
#         label='%10.1e | %6.2e' % (labeling[i][0],labeling[i][1]))

#     l = plt.legend(prop={'size':5}, title='grainsize | abundance', borderpad=0.5)
#     plt.setp(l.get_title(), fontsize='5')
#     l.get_frame().set_alpha(0.5)
#     #plt.xlim([0.03,1e4])
#     #plt.ylim([1,2e5])
#     plt.xscale('log')
#     plt.yscale('log')
#     plt.xlabel('$\lambda$ [$\mu$m]')
#     plt.ylabel('opacity [cm$^2$ g$^{-1}$]')
#     plt.title('Opacities')

    # plot kappas without labels
    kappas_name = 'kappas.dat'
    print 'kappas name'
    print kappas_name
    #select_list = [0]
    #select_list = [0,4,8,12,16,19]
    # this is a setup with 180 particles, e.g. for input from one of Til's codes.
    # each particles gets a color, and you can select which particles should be plotted using the select list
    select_list = [0, 50, 100, 150, 178]#, 4, 8, 12, 16, 19]
    print 'da3'
    color_list = ['b', 'g', 'r', 'c', 'm', 'k', 'y', 'b', 'g', 'r', 'c', 'm', 'k', 'y',
                  'b', 'g', 'r', 'c', 'm', 'k', 'y', 'b', 'g', 'r', 'c', 'm', 'k', 'y',
                  'b', 'g', 'r', 'c', 'm', 'k', 'y', 'b', 'g', 'r', 'c', 'm', 'k', 'y',
                  'b', 'g', 'r', 'c', 'm', 'k', 'y', 'b', 'g', 'r', 'c', 'm', 'k', 'y',
                  'b', 'g', 'r', 'c', 'm', 'k', 'y', 'b', 'g', 'r', 'c', 'm', 'k', 'y',
                  'b', 'g', 'r', 'c', 'm', 'k', 'y', 'b', 'g', 'r', 'c', 'm', 'k', 'y',
                  'b', 'g', 'r', 'c', 'm', 'k', 'y', 'b', 'g', 'r', 'c', 'm', 'k', 'y',
                  'b', 'g', 'r', 'c', 'm', 'k', 'y', 'b', 'g', 'r', 'c', 'm', 'k', 'y',
                  'b', 'g', 'r', 'c', 'm', 'k', 'y', 'b', 'g', 'r', 'c', 'm', 'k', 'y',
                  'b', 'g', 'r', 'c', 'm', 'k', 'y', 'b', 'g', 'r', 'c', 'm', 'k', 'y',
                  'b', 'g', 'r', 'c', 'm', 'k', 'y', 'b', 'g', 'r', 'c', 'm', 'k', 'y',
                  'b', 'g', 'r', 'c', 'm', 'k', 'y', 'b', 'g', 'r', 'c', 'm', 'k', 'y',
                  'b', 'g', 'r', 'c', 'm', 'k', 'y', 'b', 'g', 'r', 'c', 'm', 'k', 'y',
                  'b', 'g', 'r', 'c', 'm', 'k', 'y', 'b', 'g', 'r', 'c', 'm', 'k', 'y',
                  'b', 'g', 'r', 'c', 'm', 'k', 'y', 'b', 'g', 'r', 'c', 'm', 'k', 'y',
                  'b', 'g', 'r', 'c', 'm', 'k', 'y', 'b', 'g', 'r', 'c', 'm', 'k', 'y',
                  'b', 'g', 'r', 'c', 'm', 'k', 'y', 'b', 'g', 'r', 'c', 'm', 'k', 'y',
                  'b', 'g', 'r', 'c', 'm', 'k', 'y', 'b', 'g', 'r', 'c', 'm', 'k', 'y',
                  'b', 'g', 'r', 'c', 'm', 'k', 'y', 'b', 'g', 'r', 'c', 'm', 'k', 'y',
                  'b', 'g', 'r', 'c', 'm', 'k', 'y', 'b', 'g', 'r', 'c', 'm', 'k', 'y',
                  'b', 'g', 'r', 'c', 'm', 'k', 'y', 'b', 'g', 'r', 'c', 'm', 'k', 'y',
                  'b', 'g', 'r', 'c', 'm', 'k', 'y', 'b', 'g', 'r', 'c', 'm', 'k', 'y',
                  'b', 'g', 'r', 'c', 'm', 'k', 'y', 'b', 'g', 'r', 'c', 'm', 'k', 'y',
                  'b', 'g', 'r', 'c', 'm', 'k', 'y', 'b', 'g', 'r', 'c', 'm', 'k', 'y',
                  'b', 'g', 'r', 'c', 'm', 'k', 'y', 'b', 'g', 'r', 'c', 'm', 'k', 'y',
                  'b', 'g', 'r', 'c', 'm', 'k', 'y', 'b', 'g', 'r', 'c', 'm', 'k', 'y',
                  'b', 'g', 'r', 'c', 'm', 'k', 'y', 'b', 'g', 'r', 'c', 'm', 'k', 'y',
                  'b', 'g', 'r', 'c', 'm', 'k', 'y', 'b', 'g', 'r', 'c', 'm', 'k', 'y',
                  'b', 'g', 'r', 'c', 'm', 'k', 'y', 'b', 'g', 'r', 'c', 'm', 'k', 'y',
                  'b', 'g', 'r', 'c', 'm', 'k', 'y']
    #print 'da2'
    ax2 = multi_fig.add_subplot(4, 2, 2)
    #print ' da'
    wave, kappas_abs, kappas_scat = mlr.kappa_data_nolabel(path, kappas_name)
    #wave, kappas_abs, kappas_scat, labeling = mlr.selected_kappa_data(path, kappas_name)
    
    kappas_abs = np.transpose(kappas_abs)
    kappas_scat = np.transpose(kappas_scat)

#    print (wave.shape)
#    print (kappas_abs[0].shape)
    # the kappa for absorpiton and for scattering are both plotted
    j=0
    for i in xrange(len(kappas_abs)):#plot all grains
    #for i in select_list:
        ax2.plot(wave, kappas_abs[i], c= color_list[j])#, label='%10.1e | %6.1e' % (labeling[i][0],labeling[i][1]))
        j = j+1
    j = 0
    for i in xrange(len(kappas_abs)):#plot all grains
    #for i in select_list:
        ax2.plot(wave, kappas_scat[i], '--', c=color_list[j])
        j = j+1
        
#    l = plt.legend(prop={'size':5}, title='grainsize | abundance', borderpad=0.3)
#    plt.setp(l.get_title(), fontsize='4')
    #l.get_frame().set_alpha(0.5)
    #plt.xlim([0.03,1e4])
    plt.ylim([1e-3,4e5])
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('$\lambda$ [$\mu$m]')
    plt.ylabel('opacity [cm$^2$ g$^{-1}$]')
    plt.title('Opacities')

    # plot temperatur and density with tau
    denstemp_name = 'denstemp' + timestep + '.fits.gz'
    height_name = 'height' + timestep + '.dat'
    heightR_name = 'heightR' + timestep + '.dat'
    number_tau1 = 2 # assume that extra tau1 surfaces have been calculated

    radius, z_over_r, dens, temp = mlr.denstemp_data(path, denstemp_name)
    h_radius, tau_h = mlr.height_data(path, height_name)
    hR_radius, tau_hR = mlr.heightR_data(path, heightR_name)

    # comment in to add additional tau=1 surface
    #add addition tau=1 surfaces
    #hR_radius2, tau_hR2, hR_rad_list, tau_hR_list = mlr.multi_heightR_data(
    #     path, heightR_name, number_tau1)
    #h_radius2, tau_h2= mlr.multi_height_data(path, height_name, number_tau1)

    #print(radius.shape)
    #print(z_over_r.shape)
    #print(deyns.shape)

    ax3 = multi_fig.add_subplot(4, 2, 3)
    #im = ax3.pcolor(radius, z_over_r, dens, norm=LogNorm(vmin=1.e-34))
    im = ax3.pcolor(radius, z_over_r, dens, norm=LogNorm(vmin=1.e-34, vmax=1.e-8), rasterized=True)
    plt.xscale('log')
    plt.yscale('log')
#    plt.xlim([1e-1,5.e2])#DIANA27
#    plt.ylim([4.e-3,3e1])#DIANA27
    plt.xlim([1e-2,2e2])#Paola
    plt.ylim([5e-4,1e2])#Paola
#    plt.xlim([1.e-1,2e2])#Kama
#    plt.ylim([2.e-3,3.e1])#Kama
    plt.colorbar(im, aspect=10, label='density [g cm$^{-3}$]')
    plt.plot(h_radius, tau_h, label='tau=1 for $\lambda$=0.55 $\mathrm{\mu}$m')
    #plt.plot(hR_radius, tau_hR, label='tauR=1 for $\lambda$=1.67 $\mathrm{\mu}$m')
    #plt.plot(h_radius2, tau_h2[1], c='b', label='tau=1 for $\lambda$=1.67 $\mathregular{\mu}$m')
     #plt.plot(hR_rad_list[0], tau_hR_list[0], ':', c='b')
#    plt.plot(h_radius2, tau_h2[2], c='g', label='tau=1 for $\lambda$=1.6 $\mathregular{\mu}$m')
#    plt.plot(hR_rad_list[1], tau_hR_list[1], ':', c='g')
#    plt.plot(h_radius2, tau_h2[3], c='c', label='tau=1 for $\lambda$=850 $\mathregular{\mu}$m')
#    plt.plot(hR_rad_list[2], tau_hR_list[2], ':', c='c')
#    plt.plot([2e-1],[2e-3], ':', c='k', label='Rtau1, same $\lambda$')

    l = plt.legend(prop = {'size':5}, borderpad=0.5)
    l.get_frame().set_alpha(0.5)
    plt.xlabel('r [AU]')
    plt.ylabel('z/r')
    plt.title('Density')


    #plot temperature
    ax4 = multi_fig.add_subplot(4, 2, 4)
    #im = ax4.pcolor(radius, z_over_r, temp, norm=LogNorm())
    im = ax4.pcolor(radius, z_over_r, temp, norm=LogNorm(vmin=1., vmax=2.e3), rasterized=True)
    plt.xscale('log')
    plt.yscale('log')
    #plt.xlim([1e-1,5.e2])#DIANA27
    #plt.ylim([4.e-3,3e1])#DIANA27
#    plt.xlim([1.e-2,2e2])#Kama
#    plt.ylim([2.e-3,1.1])#Kama
    plt.xlim([1e-2,2e2])#Paola
    plt.ylim([5e-4,1e2])#Paola
    plt.colorbar(im, aspect=10, label='temperature [K]')
    plt.plot(h_radius, tau_h, c='k', label='tau=1 for $\lambda_{\mathregular{star}}$')
    plt.plot(hR_radius, tau_hR, ':', c='k', label='tauR=1 for $\lambda_{\mathregular{star}}$')
    l = plt.legend(prop={'size':5}, borderpad=0.5)
    l.get_frame().set_alpha(0.5)
    plt.xlabel('r [AU]')
    plt.ylabel('z/r')
    plt.title('Temperature')


    #plot surfacedensity
    surfacedens_name = 'surfacedens' + timestep + '.dat'
    ax5 = multi_fig.add_subplot(4, 2, 5)

    radius, surfacedens = mlr.surfacedens_data(path, surfacedens_name)

    ax5.plot(radius, surfacedens)
    plt.xscale('log')
    plt.yscale('log')
    #plt.xlim([1e-1,5.e2])#DIANA27
    plt.xlim([1e-2,2e2])
    plt.ylim(bottom=1e-6)
    plt.xlabel('r [AU]')
    plt.ylabel('surfacedensity [g cm$^{-2}$]')
    plt.title('Surfacedensity')

    #plot scalheights and Q
    scaleheight_name = 'scaleheight' + timestep + '.dat'
    ax6 = multi_fig.add_subplot(4, 2, 6)

    radius, h_gas, h_dust, Q_t = mlr.scaleheights_Q_data(path, scaleheight_name)
    plot_Q_t = Q_t[np.where(Q_t < 15)]
    plot_Q_t_r = radius[np.where(Q_t < 15)]

    ax6.plot(radius, h_gas, label='h$_{\mathrm{gas}}$')
    ax6.plot(radius, h_dust, label='h$_{\mathrm{dust}}$')
    ax6.plot(plot_Q_t_r, plot_Q_t, '+', label='Toomre Q < 15')
    l = plt.legend(prop={'size':5}, borderpad=0.5, loc=2)
    plt.setp(l.get_title(), fontsize='5')
    l.get_frame().set_alpha(0.5)
    plt.xscale('log')
    plt.yscale('log')
    #plt.ylim(1e-3,1e2)
    #plt.xlim([1e-1,5.e2])#DIANA27
    plt.xlim([1e-2,1e3])#Paola
    #plt.xlim([1e-2, 2e2])
    plt.xlabel('r [AU]')
    plt.ylabel('z/r')
    plt.title('Scaleheights')

    #plot spectrum for different inclinations or Visibilities
    ax7 = multi_fig.add_subplot(4, 2, 7)

    if (basevis != ''):
        
        #wave_index = 58 #corresponds to 1.68 my, line number-1 in MCSpec file
        #wave_index = 133 #Kama
        #wave_index = 78
        basevis_name = basevis
        number_wave_x_baseangle = 1 #plot all baselines for the first angle
        base, full_disk, vis, phase = mlr.basevis_pionier(path, basevis_name, number_wave_x_baseangle)
        

        # comment out add the fraction of the stellar flux and the fractoin of scattered light as a line to the visibility plot
        #all_angles = [10.5, 18.2, 23.6, 28.0, 31.8, 35.2, 38.4, 41.4, 44.2, 46.9,
        #              49.5, 51.9, 54.3, 56.6, 58.9, 61.1, 63.3, 65.4, 67.5, 69.5,
        #              71.5, 73.5, 75.5, 77.5, 79.4, 81.4, 83.3, 85.2, 87.1, 89.0]

        #angle = 44.2

        #MCSpec_wave, MCSpec_flux_dens, MCSpec_flux = mlr.MCSpec_data(path, angle)
        #MCSpec_wave2, MCSpec_flux_star, MCSpec_flux_scatter = mlr.MCSpec_star_scatter(path, angle)

        #Rel_flux_star = MCSpec_flux_star[wave_index]/MCSpec_flux[wave_index]
        #Rel_flux_scatter = MCSpec_flux_scatter[wave_index]/MCSpec_flux[wave_index]

        ##comment in to add data to the visibility plot
        #observation = 'hd100453'
        #observation = '' 
        #if observation == 'hd100453':
        #    base_obs, vis2_obs = ifq.qread_infits('/home/klarmann/HD100453_GRAVITY/pionier_data/')
        #    #print np.shape(base_obs)
        #    #print np.shape(vis2_obs)
        #label_list = [' 2.2 $\mathrm{\mu}$m', '10.7 $\mathrm{\mu}$m']
        #lambda_list = [2.2, 10.7]
        # m_list = 4 * ['+', 'o', '^']
        # ls_list = ['-', '--', ':']
        c_list = 6 * ['k', 'b']
        #label_list = ['1.60 m$\mu$ 10$\degree$', '1.68 m$\mu$ 0$\degree$',
        #          '1.76 m$\mu$ 0$\degree$'] # Kama
        ls_list = ['-', '--', ':']
        m_list = 3 * ['x', 'o', '^', '+']
        #lambda_list = [1.6, 1.68, 1.76]
        #c_list = ['k', 'b', 'y']
        #if observation == 'hd100453':
            #for i in xrange(22): #number of obseration files                
                #ax7.plot(base_obs[i], vis2_obs[i][0,:],'o', c='beige', markeredgecolor='beige')
                #ax7.plot(base_obs[i], vis2_obs[i][2,:], 'o', c='snow', markeredgecolor='snow')
                #ax7.plot(base_obs[i]/1.68, vis2_obs[i][1,:], 'o', c='yellow',  markersize=2.9)
                #plt.hlines(0.36, 0, 100, colors=u'r', linestyles=':')
                #plt.vlines(5, 0, 1, colors=u'k', linestyles='--')
        for i in [0]:
            print i
            rel_base = base # divide by baseline for baseline in Mlambda
            vis2 = []
            for j in xrange(len(vis[:,i])):
                vis2.append(vis[j,i] * vis[j,i])
            print len(vis2)
            ax7.plot(rel_base, vis2, ls=ls_list[i])

        # comment in for rel flux from star and scattered light
        #plt.hlines(Rel_flux_star**2, base[0], base[-1],
        #           colors=u'r', linestyles=u'solid', label=u'star')
        #print Rel_flux_star
        #plt.hlines(1.-Rel_flux_scatter**2, base[0], base[-1],
        #           colors=u'g', linestyles=u'solid', label=u'scatter')
        #print Rel_flux_scatter
               
        ax7.set_xlim([0,130])
        ax7.set_ylim([0, 1.])
        ax7.set_xlabel('baseline [m]')
        ax7.set_ylabel('v2')
        #ax7.legend(prop = {'size':6}, borderpad=0.5)
        ax7.set_title('Visibility squared ')

    # if there is no basevis file, plot sed for different angles    
    else:
        all_angles = [10.5, 18.2, 23.6, 28.0, 31.8, 35.2, 38.4, 41.4, 44.2, 46.9,
                      49.5, 51.9, 54.3, 56.6, 58.9, 61.1, 63.3, 65.4, 67.5, 69.5,
                      71.5, 73.5, 75.5, 77.5, 79.4, 81.4, 83.3, 85.2, 87.1, 89.0]

        angle_index = [0, 5, 10, 15, 20, 25, 29]

        angles = [all_angles[int(i)]for i in angle_index]
    #    print(angles)

        for angle in angles:
            MCSpec_wave, MCSpec_flux_dens, MCSpec_flux = mlr.MCSpec_data(path, angle)
            ax7.plot(MCSpec_wave, MCSpec_flux, label=str(angle)+' deg')

        l = plt.legend(prop={'size':5}, title='inclination', borderpad=0.5, ncol=3, loc=8)
        plt.setp(l.get_title(), fontsize='5')
        l.get_frame().set_alpha(0.5)
        plt.xscale('log')
        plt.yscale('log')
        plt.xlabel('$\lambda$ [$ \mathregular{\mu}$m]')
        plt.ylabel('$\lambda F_\lambda$ [erg cm$^{-2}$s$^{-1}$]')
        plt.ylim([1e-12,5e-8])
        plt.title('Inclinations')

#plot temperature and density linear zoom in 

    radius, z_over_r, dens, temp = mlr.denstemp_data(path, denstemp_name, unit='AU')
    h_radius, tau_h = mlr.height_data(path, height_name, unit='AU')
    hR_radius, tau_hR = mlr.heightR_data(path, heightR_name, unit='AU')

    ax8 = multi_fig.add_subplot(4, 2, 8)
    #im = ax8.pcolor(radius, z_over_r, temp, norm=LogNorm())
    im = ax8.pcolor(radius, z_over_r, temp, norm=LogNorm(vmin=1., vmax=2.e3), rasterized=True)
    #plt.xscale('log')
    #plt.yscale('log')
    plt.plot(h_radius, tau_h, c='k', label='tau=1 for $\lambda_{\mathregular{star}}$')
    plt.plot(hR_radius, tau_hR, ':', c='k', label='tauR=1 for $\lambda_{\mathregular{star}}$')
    plt.xlim([0.,5.])
    plt.ylim([0.,1.])
#    plt.xlim([5.5e-1,5e2])
#    plt.ylim([2.e-3,1.1])
    plt.colorbar(im, aspect=10, label='temperature [K]')
#    plt.plot(h_radius, tau_h, c='k', label='tau=1 for $\lambda_{\mathregular{star}}$')
#    plt.plot(hR_radius, tau_hR, ':', c='k', label='tauR=1 for $\lambda_{\mathregular{star}}$')
    l = plt.legend(prop={'size':4}, borderpad=0.4)
    l.get_frame().set_alpha(0.5)
    plt.xlabel('r [AU]')
    plt.ylabel('z [AU]')
    plt.title('Temperature')


    a = plt.axes([.65, .15, .095, .045])
    #n, bins, patches = plt.hist(s, 400, normed=1)
    im2 = a.pcolor(radius, z_over_r, dens, norm=LogNorm(vmin=1.e-34, vmax=1.e-8), rasterized=True)
    plt.plot(h_radius, tau_h, c='k', label='tau=1 for $\lambda_{\mathregular{star}}$')
    plt.plot(hR_radius, tau_hR, ':', c='k', label='tauR=1 for $\lambda_{\mathregular{star}}$')
    cb = plt.colorbar(im2, aspect=10, ticks=[1.e-34, 1.e-27, 1.e-21, 1.e-14, 1.e-8],)#, label='density [g cm$^{-3}$]')
     # grab the Colorbar instance
    for t in cb.ax.get_yticklabels():
        t.set_fontsize(5)
    plt.xlim([0.,5.])
    plt.ylim([0.,1.])
    plt.xticks(np.arange(6))
    plt.tick_params(axis='both', which='major', labelsize=5)
    plt.title('Density', size=7)



    print 'before tight'
    #multi_fig.suptitle('selfconsitent setteling')
#    plt.text(0.5, 0.03, 'hotter disks', transform=multi_fig.transFigure, horizontalalignment='center')
    multi_fig.tight_layout()
    print 'first save'
    #save to 'plot' folder inside the Output folder
    plt.savefig(plotpath_multi+model_name+"_"+timestep+'multi.pdf',
                format='pdf', dpi=200, bbox_inches='tight')
    print 'second save'
    # save to common 'plots' folder, located in the same folder where all the model folders are located
    plt.savefig(common_path[0:-1] + '/plots/' +  model_name + '_multi.pdf',
                format='pdf', dpi=200, bbox_inches='tight')  

    plt.close('all')
    return

