import numpy as np
import pygad
import datetime as dt
from spacepy import pycdf
import bisect
import scipy.constants as const
import pandas as pd
from scipy.interpolate import griddata as gd
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import ticker
from matplotlib.ticker import LogFormatter


def plot_original_vdf(dir_fig, filename, azimuth, elevation, velocity, vdf):
    fig, axs = plt.subplots(3, 3, subplot_kw={'projection': 'polar'}, figsize=(9, 8))
    plt.subplots_adjust(wspace=0.05, hspace=0.3, left=0.01, right=0.9, top=0.9, bottom=0.05)
    axs = axs.flatten()
    # ax = axs[1 - 1]
    for iax, ax in enumerate(axs):
        if iax < 9:
            img_arr = vdf[:, iax, :]
            img_arr_plot = np.ma.array(img_arr, mask=img_arr == 0)
            x, y = np.meshgrid(np.radians(azimuth), velocity)
            c = np.ones_like(x)
            ax.grid(False)
            ax.pcolormesh(x, y, c, facecolor='none', edgecolor='k', linewidth=0.01, linestyle=':')
            # print(np.min(img_arr_plot),np.max(img_arr_plot))
            img = ax.pcolormesh(np.radians(azimuth), velocity, img_arr_plot.T, cmap='jet', shading='auto',
                                norm=mpl.colors.LogNorm(
                                    vmin=np.min(np.ma.array(vdf, mask=np.logical_or(vdf == 0, velocity < 250))),
                                    vmax=np.max(np.ma.array(vdf, mask=np.logical_or(vdf == 0, velocity < 250)))))

            ax.set_ylim([250, 600])
            ax.set_title('elevation=' + '{:.1f}'.format(elevation[iax]) + r'$^{\circ}$')
            ax.xaxis.set_label_coords(1.0, -0.05)
            ax.set_xlabel('(km/s)')
            ax.set_rorigin(0)
            ax.xaxis.set_major_locator(ticker.MultipleLocator(np.radians(10)))
            ax.yaxis.set_major_locator(ticker.MultipleLocator(100))
            ax.set_xlim([np.min(np.radians(azimuth)), np.max(np.radians(azimuth))])
            if iax == 2:
                fig.colorbar(img, ax=ax, label=r'VDF $(cm^{-3}\cdot s^6)$')
    fig.suptitle(str(epoch))
    plt.savefig(dir_fig + filename, dpi=300)
    plt.close()


def plot_rtn_mfa_vdf(dir_fig, filename, grid_vdf_3d_rtn, grid_vdf_3d_mfa, max_ind_rtn, max_ind_mfa,
                     vr_lst, vt_lst, vn_lst, vpar_lst, vper1_lst, vper2_lst):
    fig, axs = plt.subplots(2, 3, figsize=(15, 10))
    plt.subplots_adjust(wspace=0.25, hspace=0.3, left=0.08, right=0.90, top=0.9, bottom=0.1)
    axs = axs.flatten()
    # set vmin and vmax
    grid_vdf_3d_plot = np.copy(grid_vdf_3d_rtn)
    grid_vdf_3d_plot[np.where(grid_vdf_3d_plot == 0.0)] = np.nan
    img_plot_1 = grid_vdf_3d_plot[:, :, max_ind_rtn[2]].T
    img_plot_2 = grid_vdf_3d_plot[:, max_ind_rtn[1], :].T
    img_plot_3 = grid_vdf_3d_plot[max_ind_rtn[0], :, :].T
    grid_vdf_3d_plot = np.copy(grid_vdf_3d_mfa)
    grid_vdf_3d_plot[np.where(grid_vdf_3d_plot == 0.0)] = np.nan
    img_plot_4 = grid_vdf_3d_plot[:, :, max_ind_mfa[2]].T
    img_plot_5 = grid_vdf_3d_plot[:, max_ind_mfa[1], :].T
    img_plot_6 = grid_vdf_3d_plot[max_ind_mfa[0], :, :].T
    val_rtn_min = np.nanpercentile([img_plot_1, img_plot_2, img_plot_3], 0)
    val_rtn_max = np.nanpercentile([img_plot_1, img_plot_2, img_plot_3], 100)
    val_mfa_min = np.nanpercentile([img_plot_4, img_plot_5, img_plot_6], 0)
    val_mfa_max = np.nanpercentile([img_plot_4, img_plot_5, img_plot_6], 100)
    levels_rtn = np.logspace(np.log10(val_rtn_min), np.log10(val_rtn_max), 32)
    levels_mfa = np.logspace(np.log10(val_mfa_min), np.log10(val_mfa_max), 32)
    vr_max, vt_max, vn_max = vr_lst[max_ind_rtn[0]], vt_lst[max_ind_rtn[1]], vn_lst[max_ind_rtn[2]]
    vr_beg, vt_beg, vn_beg = vr_max - epara[0] * 100, vt_max - epara[1] * 100, vn_max - epara[2] * 100

    ax = axs[1 - 1]
    img_plot = img_plot_1
    xlim = [np.min(vr_lst), np.max(vr_lst)]
    ylim = [np.min(vt_lst), np.max(vt_lst)]
    img_arr = np.ma.masked_array(img_plot, mask=np.isnan(img_plot))
    img = ax.contourf(vr_lst, vt_lst, img_arr, levels=levels_rtn,
                      norm=mpl.colors.LogNorm(vmin=val_rtn_min, vmax=val_rtn_max), cmap='jet')
    ax.arrow(vr_beg, vt_beg, 200 * epara[0], 200 * epara[1], width=1, head_width=5)
    ax.set_xlim(xlim[0], xlim[1])
    ax.set_ylim(ylim[0], ylim[1])
    ax.set_xlabel(r'$V_{R} (km/s)$', fontsize=15)
    ax.set_ylabel(r'$V_{T} (km/s)$', fontsize=15)
    ax.set_aspect('equal')
    ax.grid(True)

    ax = axs[2 - 1]
    img_plot = img_plot_2
    xlim = [np.min(vr_lst), np.max(vr_lst)]
    ylim = [np.min(vn_lst), np.max(vn_lst)]
    img_arr = np.ma.masked_array(img_plot, mask=np.isnan(img_plot))
    img = ax.contourf(vr_lst, vn_lst, img_arr, levels=levels_rtn,
                      norm=mpl.colors.LogNorm(vmin=val_rtn_min, vmax=val_rtn_max), cmap='jet')
    ax.arrow(vr_beg, vn_beg, 200 * epara[0], 200 * epara[2], width=1, head_width=5)
    ax.set_xlim(xlim[0], xlim[1])
    ax.set_ylim(ylim[0], ylim[1])
    ax.set_xlabel(r'$V_{R} (km/s)$', fontsize=15)
    ax.set_ylabel(r'$V_{N} (km/s)$', fontsize=15)
    ax.set_aspect('equal')
    ax.grid(True)

    ax = axs[3 - 1]
    img_plot = img_plot_3
    xlim = [np.min(vt_lst), np.max(vt_lst)]
    ylim = [np.min(vn_lst), np.max(vn_lst)]
    img_arr = np.ma.masked_array(img_plot, mask=np.isnan(img_plot))
    img = ax.contourf(vt_lst, vn_lst, img_arr, levels=levels_rtn,
                      norm=mpl.colors.LogNorm(vmin=val_rtn_min, vmax=val_rtn_max), cmap='jet')
    ax.arrow(vt_beg, vn_beg, 200 * epara[1], 200 * epara[2], width=1, head_width=5)
    ax.set_xlim(xlim[0], xlim[1])
    ax.set_ylim(ylim[0], ylim[1])
    ax.set_xlabel(r'$V_{T} (km/s)$', fontsize=15)
    ax.set_ylabel(r'$V_{N} (km/s)$', fontsize=15)
    ax.set_aspect('equal')
    ax.grid(True)
    cb_ax = fig.add_axes([.92, .55, .01, .35])
    fig.colorbar(img, orientation='vertical', cax=cb_ax, label=r'VDF ($m^{-3}\cdot s^6$)')

    ax = axs[4 - 1]
    img_plot = img_plot_4
    xlim = [np.min(vpar_lst), np.max(vpar_lst)]
    ylim = [np.min(vper1_lst), np.max(vper1_lst)]
    img_arr = np.ma.masked_array(img_plot, mask=np.isnan(img_plot))
    img = ax.contourf(vpar_lst, vper1_lst, img_arr, levels=levels_mfa,
                      norm=mpl.colors.LogNorm(vmin=val_mfa_min, vmax=val_mfa_max), cmap='jet')
    ax.set_xlim(xlim[0], xlim[1])
    ax.set_ylim(ylim[0], ylim[1])
    ax.set_xlabel(r'$V_{||} (km/s)$', fontsize=15)
    ax.set_ylabel(r'$V_{\perp,1} (km/s)$', fontsize=15)
    ax.set_aspect('equal')
    ax.grid(True)

    ax = axs[5 - 1]
    img_plot = img_plot_5
    xlim = [np.min(vpar_lst), np.max(vpar_lst)]
    ylim = [np.min(vper2_lst), np.max(vper2_lst)]
    img_arr = np.ma.masked_array(img_plot, mask=np.isnan(img_plot))
    img = ax.contourf(vpar_lst, vper2_lst, img_arr, levels=levels_mfa,
                      norm=mpl.colors.LogNorm(vmin=val_mfa_min, vmax=val_mfa_max), cmap='jet')
    ax.set_xlim(xlim[0], xlim[1])
    ax.set_ylim(ylim[0], ylim[1])
    ax.set_xlabel(r'$V_{||} (km/s)$', fontsize=15)
    ax.set_ylabel(r'$V_{\perp,2} (km/s)$', fontsize=15)
    ax.set_aspect('equal')
    ax.grid(True)

    ax = axs[6 - 1]
    img_plot = img_plot_6
    xlim = [np.min(vper1_lst), np.max(vper1_lst)]
    ylim = [np.min(vper2_lst), np.max(vper2_lst)]
    img_arr = np.ma.masked_array(img_plot, mask=np.isnan(img_plot))
    img = ax.contourf(vper1_lst, vper2_lst, img_arr, levels=levels_mfa,
                      norm=mpl.colors.LogNorm(vmin=val_mfa_min, vmax=val_mfa_max), cmap='jet')
    ax.set_xlim(xlim[0], xlim[1])
    ax.set_ylim(ylim[0], ylim[1])
    ax.set_xlabel(r'$V_{\perp,1} (km/s)$', fontsize=15)
    ax.set_ylabel(r'$V_{\perp,2} (km/s)$', fontsize=15)
    ax.set_aspect('equal')
    ax.grid(True)
    cb_ax = fig.add_axes([.92, .1, .01, .35])
    fig.colorbar(img, orientation='vertical', cax=cb_ax, label=r'VDF ($m^{-3}\cdot s^6$)')
    # plt.show()
    plt.savefig(dir_fig + filename, dpi=300)
    plt.close()


def fitness_func(solution, solution_idx):
    output = 1.e-3 * solution[0] / (solution[2] ** 2 * solution[3] * np.pi ** 1.5) * \
             np.exp(-(vpara - solution[6]) ** 2 / solution[3] ** 2) * \
             np.exp(-(vperp1 - solution[7]) ** 2 / solution[2] ** 2) * \
             np.exp(-(vperp2 - solution[8]) ** 2 / solution[2] ** 2) + \
             1.e-3 * solution[1] / (solution[4] ** 2 * solution[5] * np.pi ** 1.5) * \
             np.exp(-(vpara - solution[9]) ** 2 / solution[5] ** 2) * \
             np.exp(-(vperp1 - solution[10]) ** 2 / solution[4] ** 2) * \
             np.exp(-(vperp2 - solution[11]) ** 2 / solution[4] ** 2)
    # fitness = 1.e-5 / np.nansum(np.sqrt(np.abs((output - desired_output)) ** 2.0))
    fitness = -np.nansum(np.sqrt(np.abs((output - desired_output)) ** 2.0))
    # fitness = -np.nansum(np.abs((np.log(output) - np.log(desired_output))) ** 2.0)
    return fitness


def on_generation(ga_instance):
    global last_fitness
    print("Generation = {generation}".format(generation=ga_instance.generations_completed))
    print("Fitness    = {fitness}".format(
        fitness=ga_instance.best_solution(pop_fitness=ga_instance.last_generation_fitness)[1]))
    print("Change     = {change}".format(
        change=ga_instance.best_solution(pop_fitness=ga_instance.last_generation_fitness)[1] - last_fitness))
    last_fitness = ga_instance.best_solution(pop_fitness=ga_instance.last_generation_fitness)[1]


# set calculation time
year, mon, day = 2020, 8, 2
_year, _mon, _day = str(year), str(mon), str(day)
beg_time = dt.datetime(year, mon, day, 11, 38, 34, 000000)
end_time = dt.datetime(year, mon, day, 11, 38, 35, 000000)
_date = beg_time.strftime('%Y%m%d')
_beg_t, _end_t = beg_time.strftime('%H%M%S.%f'), end_time.strftime('%H%M%S.%f')
dir_cdf = '/Users/psr/work/pku/nfa_beam/Data/'
# read proton moments cdf file
pas_file = 'solo_l2_swa-pas-vdf_20200802_v02.cdf'
pas_cdf = pycdf.CDF(dir_cdf + pas_file)
print(pas_cdf)
start_ind = bisect.bisect_left(pas_cdf['Epoch'], beg_time)
stop_ind = bisect.bisect_right(pas_cdf['Epoch'], end_time)
epoch_lst = pas_cdf['Epoch'][start_ind:stop_ind]
azimuth, elevation, energy = np.array(pas_cdf['Azimuth']), np.array(pas_cdf['Elevation']), np.array(pas_cdf['Energy'])
velocity = np.sqrt(2 * energy * const.e / const.m_p) * 1.e-3
vdf_lst = pas_cdf['vdf'][start_ind:stop_ind, :, :, :]
pas_to_rtn_lst = pas_cdf['PAS_to_RTN'][start_ind:stop_ind, :, :]
# read magnetic field cdf file
mag_file = 'solo_l2_mag-rtn-normal_20200802_v02.cdf'
mag_cdf = pycdf.CDF(dir_cdf + mag_file)
print(mag_cdf)
start_ind = bisect.bisect_left(mag_cdf['EPOCH'], beg_time)
stop_ind = bisect.bisect_right(mag_cdf['EPOCH'], end_time)
mag_epoch = mag_cdf['EPOCH'][start_ind:stop_ind]
br_lst = mag_cdf['B_RTN'][start_ind:stop_ind, 0]
bt_lst = mag_cdf['B_RTN'][start_ind:stop_ind, 1]
bn_lst = mag_cdf['B_RTN'][start_ind:stop_ind, 2]
# interpolation
julday_lst = np.array(pd.to_datetime(epoch_lst).to_julian_date())
julday_mag_lst = np.array(pd.to_datetime(mag_epoch).to_julian_date())
br_interp = np.interp(julday_lst, julday_mag_lst, br_lst)
bt_interp = np.interp(julday_lst, julday_mag_lst, bt_lst)
bn_interp = np.interp(julday_lst, julday_mag_lst, bn_lst)

ele_arr = np.tile(elevation[None, :, None], (11, 1, 96))
azi_arr = np.tile(azimuth[:, None, None], (1, 9, 96))
energy_arr = np.tile(energy[None, None, :], (11, 9, 1))
vel_arr = np.sqrt(2 * energy_arr * const.e / const.m_p) * 1.e-3
# whether to plot original vdf or conduct genetic algorithm fitting
is_plot_vdf_in_sph_mesh = 0
is_plot_original_vdf = 0
for i_time, epoch in enumerate(epoch_lst):
    vdf = vdf_lst[i_time, :, :, :]
    # plot original vdf (this plot must be ahead of the following)
    # ===================================================
    if is_plot_vdf_in_sph_mesh == 1:
        dir_fig = '/Users/psr/work/pku/nfa_beam/Figures/original_vdf/' + epoch.strftime('%Y-%m-%d') + '/'
        filename = 'mesh_vdf_pas_frame_' + epoch.strftime('%H%M%S.%f')[:-3] + '.png'
        plot_original_vdf(dir_fig, filename, azimuth, elevation, velocity, vdf)
    # ===================================================
    # set vdfs at velocities less than 300 km/s to 0.0 which are not reliable
    if dt.datetime(2020, 8, 2, 11, 38, 34) < epoch < dt.datetime(2020, 8, 2, 11, 38, 35):
        vdf[np.where(vel_arr < 310)] = 0.0
    else:
        vdf[np.where(vel_arr < 300)] = 0.0
    pas_to_rtn = pas_to_rtn_lst[i_time, :, :]
    bx, by, bz = br_interp[i_time], bt_interp[i_time], bn_interp[i_time]
    vx = -np.cos(np.deg2rad(ele_arr)) * np.cos(np.deg2rad(azi_arr)) * vel_arr
    vy = -np.cos(np.deg2rad(ele_arr)) * np.sin(np.deg2rad(azi_arr)) * vel_arr
    vz = np.sin(np.deg2rad(ele_arr)) * vel_arr
    vr = vx * pas_to_rtn[0, 0] + vy * pas_to_rtn[0, 1] + vz * pas_to_rtn[0, 2]
    vt = vx * pas_to_rtn[1, 0] + vy * pas_to_rtn[1, 1] + vz * pas_to_rtn[1, 2]
    vn = vx * pas_to_rtn[2, 0] + vy * pas_to_rtn[2, 1] + vz * pas_to_rtn[2, 2]
    epara = [bx, by, bz]
    epara = epara / np.linalg.norm(epara)
    eperp2 = np.cross(epara, [1, 0, 0])
    eperp2 = eperp2 / np.linalg.norm(eperp2)
    eperp1 = np.cross(eperp2, epara)
    eperp1 = eperp1 / np.linalg.norm(eperp1)
    vpara_rtn = vr * epara[0] + vt * epara[1] + vn * epara[2]
    vperp1_rtn = vr * eperp1[0] + vt * eperp1[1] + vn * eperp1[2]
    vperp2_rtn = vr * eperp2[0] + vt * eperp2[1] + vn * eperp2[2]

    # find the maximum vdf in the original array
    # vdf[np.where(np.logical_and(np.logical_or(np.abs(vperp1_rtn) > 150, np.abs(vperp2_rtn > 150)), vdf > 1.e-9))] = 0.0
    max_ind = np.nanargmax(vdf)
    max_vdf_tmp = np.nanmax(vdf[np.logical_and(vel_arr > 300, vel_arr < 350)])
    if dt.datetime(2020, 8, 2, 11, 38, 34) < epoch < dt.datetime(2020, 8, 2, 11, 38, 35):
        max_vdf_tmp = np.nanmax(vdf[np.logical_and(vel_arr > 310, vel_arr < 350)])
    max_ind = np.ravel_multi_index(np.where(vdf == max_vdf_tmp), (11, 9, 96))
    if dt.datetime(2020, 8, 2, 11, 57, 58) < epoch < dt.datetime(2020, 8, 2, 11, 57, 59): max_ind = 3997
    if dt.datetime(2020, 8, 2, 11, 58, 18) < epoch < dt.datetime(2020, 8, 2, 11, 58, 19): max_ind = 3996
    if dt.datetime(2020, 8, 2, 11, 58, 50) < epoch < dt.datetime(2020, 8, 2, 11, 58, 51): max_ind = 4859
    # if dt.datetime(2020, 8, 2, 11, 59, 22) < epoch < dt.datetime(2020, 8, 2, 11, 59, 23): max_ind = 4861
    vpara = vpara_rtn - vpara_rtn[np.unravel_index(max_ind, (11, 9, 96))]
    vperp1 = vperp1_rtn - vperp1_rtn[np.unravel_index(max_ind, (11, 9, 96))]
    vperp2 = vperp2_rtn - vperp2_rtn[np.unravel_index(max_ind, (11, 9, 96))]

    # interpolate to 3d regular grid (rtn coordinates)
    if is_plot_original_vdf == 1:
        points = np.vstack((vr.reshape(-1), vt.reshape(-1), vn.reshape(-1))).T
        vr_lst, vt_lst, vn_lst = np.linspace(200, 600, 51), np.linspace(-200, 200, 51), np.linspace(-200, 200, 51)
        vr_arr, vt_arr, vn_arr = np.meshgrid(vr_lst, vt_lst, vn_lst, indexing='ij')
        xi = np.vstack((vr_arr.reshape(-1), vt_arr.reshape(-1), vn_arr.reshape(-1))).T
        grid_vdf_rtn = gd(points, vdf.reshape(-1), xi, method='linear')
        grid_vdf_3d_rtn = grid_vdf_rtn.reshape((51, 51, 51))
        max_ind_rtn = np.unravel_index(np.nanargmax(grid_vdf_3d_rtn), (51, 51, 51))
    # interpolate to 3d regular grid (mfa coordinates)
    points = np.vstack((vpara.reshape(-1), vperp1.reshape(-1), vperp2.reshape(-1))).T
    if is_plot_original_vdf == 1:
        vpar_lst, vper1_lst, vper2_lst = np.linspace(-350, 50, 51), np.linspace(-200, 200, 51), np.linspace(-200, 200,
                                                                                                            51)
    elif is_plot_original_vdf != 1:
        vpar_lst, vper1_lst, vper2_lst = np.linspace(-150, 50, 51), np.linspace(-100, 100, 51), np.linspace(-100, 100,
                                                                                                            51)
    vpar_arr, vper1_arr, vper2_arr = np.meshgrid(vpar_lst, vper1_lst, vper2_lst, indexing='ij')
    xi = np.vstack((vpar_arr.reshape(-1), vper1_arr.reshape(-1), vper2_arr.reshape(-1))).T
    grid_vdf_mfa = gd(points, vdf.reshape(-1), xi, method='linear', fill_value=0.0)
    grid_vdf_3d_mfa = grid_vdf_mfa.reshape((51, 51, 51))
    max_ind_mfa = np.unravel_index(np.argmax(grid_vdf_3d_mfa), (51, 51, 51))

    # # plot interpolated regular mesh
    # ===================================================
    if is_plot_original_vdf == 1:
        dir_fig = '/Users/psr/work/pku/nfa_beam/Figures/regular_int_vdf/' + epoch.strftime('%Y-%m-%d') + '/'
        filename = 'mesh_vdf_reg_rtn_frame_' + epoch.strftime('%H%M%S.%f')[:-3] + '.png'
        plot_rtn_mfa_vdf(dir_fig, filename, grid_vdf_3d_rtn, grid_vdf_3d_mfa, max_ind_rtn, max_ind_mfa,
                         vr_lst, vt_lst, vn_lst, vpar_lst, vper1_lst, vper2_lst)
    # ===================================================
    elif is_plot_original_vdf != 1:
        # set up pygad object
        # desired_output = vdf  # Function output.
        # prepare for fitting using genetic algorithm
        vdf_v2 = np.copy(vdf)
        vdf_v2[np.where(vdf < 1.e-11)] = np.nan
        desired_output = vdf_v2

        num_generations = 100  # Number of generations.
        num_parents_mating = 100  # Number of solutions to be selected as parents in the mating pool.
        sol_per_pop = 3000  # Number of solutions in the population.
        num_genes = 12  # Number of parameters

        genes_spac = [{'low': 1, 'high': 3.0}, {'low': 2.5, 'high': 8},
                      {'low':  10, 'high': 30}, {'low': 10, 'high': 20}, {'low': 30, 'high': 80}, {'low': 20, 'high': 50},
                      {'low': -10, 'high': 10}, {'low': -10, 'high': 10}, {'low': -10, 'high': 10},
                      {'low': -70, 'high': -20}, {'low': -25, 'high': 25}, {'low': -25, 'high': 25}]

        # ====== EXAMPLE: from pygad =====
        # function_inputs = [4, -2, 3.5, 5, -11, -4.7]  # Function inputs.
        # desired_output = 44  # Function output.
        #
        #
        # def fitness_func(so lution, solution_idx):
        #     output = np.sum(solution * function_inputs)
        #     fitness = 1.0 / (np.abs(output - desired_output) + 0.000001)
        #     return fitness
        #
        #
        # num_generations = 100  # Number of generations.
        # num_parents_mating = 10  # Number of solutions to be selected as parents in the mating pool.
        #
        # sol_per_pop = 20  # Number of solutions in the population.
        # num_genes = len(function_inputs)

        last_fitness = 0

        ga_instance = pygad.GA(num_generations=num_generations,
                               num_parents_mating=num_parents_mating,
                               sol_per_pop=sol_per_pop,
                               num_genes=num_genes,
                               gene_space=genes_spac,
                               fitness_func=fitness_func,
                               on_generation=on_generation,
                               parent_selection_type='sss')

        # Running the GA to optimize the parameters of the function.
        ga_instance.run()

        # Returning the details of the best solution.
        solution, solution_fitness, solution_idx = ga_instance.best_solution(ga_instance.last_generation_fitness)
        print("Parameters of the best solution : {solution}".format(solution=solution))
        print("Fitness value of the best solution = {solution_fitness}".format(solution_fitness=solution_fitness))
        print("Index of the best solution : {solution_idx}".format(solution_idx=solution_idx))

        # prediction = np.sum(np.array(function_inputs) * solution)
        # print("Predicted output based on the best solution : {prediction}".format(prediction=prediction))

        if ga_instance.best_solution_generation != -1:
            print("Best fitness value reached after {best_solution_generation} generations.".format(
                best_solution_generation=ga_instance.best_solution_generation))

        grid_vdf_3d_fit = 1.e-3 * solution[0] / (solution[2] ** 2 * solution[3] * np.pi ** 1.5) * \
                          np.exp(-(vpar_arr - solution[6]) ** 2 / solution[3] ** 2) * \
                          np.exp(-(vper1_arr - solution[7]) ** 2 / solution[2] ** 2) * \
                          np.exp(-(vper2_arr - solution[8]) ** 2 / solution[2] ** 2) + \
                          1.e-3 * solution[1] / (solution[4] ** 2 * solution[5] * np.pi ** 1.5) * \
                          np.exp(-(vpar_arr - solution[9]) ** 2 / solution[5] ** 2) * \
                          np.exp(-(vper1_arr - solution[10]) ** 2 / solution[4] ** 2) * \
                          np.exp(-(vper2_arr - solution[11]) ** 2 / solution[4] ** 2)
        grid_vdf_pc_3d_fit = 1.e-3 * solution[0] / (solution[2] ** 2 * solution[3] * np.pi ** 1.5) * \
                             np.exp(-(vpar_arr - solution[6]) ** 2 / solution[3] ** 2) * \
                             np.exp(-(vper1_arr - solution[7]) ** 2 / solution[2] ** 2) * \
                             np.exp(-(vper2_arr - solution[8]) ** 2 / solution[2] ** 2)
        grid_vdf_pb_3d_fit = 1.e-3 * solution[1] / (solution[4] ** 2 * solution[5] * np.pi ** 1.5) * \
                             np.exp(-(vpar_arr - solution[9]) ** 2 / solution[5] ** 2) * \
                             np.exp(-(vper1_arr - solution[10]) ** 2 / solution[4] ** 2) * \
                             np.exp(-(vper2_arr - solution[11]) ** 2 / solution[4] ** 2)

        # save the solutions
        dir_solution = '/Users/psr/work/pku/nfa_beam/Data/solutions/' + epoch.strftime('%Y-%m-%d') + '/'
        filename = 'solutions_' + epoch.strftime('%H%M%S.%f')[:-3] + '.txt'
        np.savetxt(dir_solution + filename, solution, fmt='%.10f')
        # save fitness value of the best solution
        dir_solution = '/Users/psr/work/pku/nfa_beam/Data/fitness/' + epoch.strftime('%Y-%m-%d') + '/'
        filename = 'fitness_value_' + epoch.strftime('%H%M%S.%f')[:-3] + '.txt'
        np.savetxt(dir_solution + filename, solution_fitness.reshape(1, -1), fmt='%.10f')

        # plot the GA instance.
        dir_fig = '/Users/psr/work/pku/nfa_beam/Figures/fitness/' + epoch.strftime('%Y-%m-%d') + '/'
        filename = 'Fitness_Function_' + epoch.strftime('%H%M%S.%f')[:-3] + '.png'
        fig = plt.figure(figsize=[6, 5])
        plt.subplots_adjust(wspace=0.05, hspace=0.3, left=0.15, right=0.9, top=0.9, bottom=0.1)
        plt.plot(ga_instance.best_solutions_fitness, linewidth=3, color='#3870FF')
        plt.title('PyGAD - Generation vs. Fitness', fontsize=14)
        plt.xlabel('Generation', fontsize=14)
        plt.ylabel('Fitness', fontsize=14)
        plt.savefig(dir_fig + filename, dpi=300)
        plt.close()

        # plot and compare fitting results
        dir_fig = '/Users/psr/work/pku/nfa_beam/Figures/compare_obs_fit_vdf/' + epoch.strftime('%Y-%m-%d') + '/'
        filename = 'obs_fit_vdf_compare_' + epoch.strftime('%H%M%S.%f')[:-3] + '.png'
        fig, axs = plt.subplots(3, 3, figsize=(15, 15))
        plt.subplots_adjust(wspace=0.25, hspace=0.3, left=0.08, right=0.90, top=0.9, bottom=0.1)
        axs = axs.flatten()
        # set vmin and vmax
        vdf_obs = np.copy(grid_vdf_3d_mfa)
        vdf_obs[np.where(vdf_obs == 0.0)] = np.nan
        img_plot_1 = vdf_obs[:, :, max_ind_mfa[2]].T
        img_plot_2 = vdf_obs[:, max_ind_mfa[1], :].T
        img_plot_3 = vdf_obs[max_ind_mfa[0], :, :].T
        vdf_fit = np.copy(grid_vdf_3d_fit)
        vdf_fit[np.where(vdf_fit == 0.0)] = np.nan
        img_plot_4 = vdf_fit[:, :, max_ind_mfa[2]].T
        img_plot_5 = vdf_fit[:, max_ind_mfa[1], :].T
        img_plot_6 = vdf_fit[max_ind_mfa[0], :, :].T
        val_min = np.nanpercentile([img_plot_1, img_plot_2, img_plot_3], 0)
        val_max = np.nanpercentile([img_plot_1, img_plot_2, img_plot_3], 100)
        levels = np.logspace(np.log10(val_min), np.log10(val_max), 32)
        vdf_pc_fit = np.copy(grid_vdf_pc_3d_fit)
        vdf_pc_fit[np.where(vdf_pc_fit == 0.0)] = np.nan
        vdf_pb_fit = np.copy(grid_vdf_pb_3d_fit)
        vdf_pb_fit[np.where(vdf_pb_fit == 0.0)] = np.nan
        lin_plot_1 = vdf_obs[:, max_ind_mfa[1], max_ind_mfa[2]]
        lin_plot_2 = vdf_obs[max_ind_mfa[0], :, max_ind_mfa[2]]
        lin_plot_3 = vdf_obs[max_ind_mfa[0], max_ind_mfa[1], :]
        lin_plot_4 = vdf_fit[:, max_ind_mfa[1], max_ind_mfa[2]]
        lin_plot_5 = vdf_fit[max_ind_mfa[0], :, max_ind_mfa[2]]
        lin_plot_6 = vdf_fit[max_ind_mfa[0], max_ind_mfa[1], :]
        lin_plot_pc_4 = vdf_pc_fit[:, max_ind_mfa[1], max_ind_mfa[2]]
        lin_plot_pc_5 = vdf_pc_fit[max_ind_mfa[0], :, max_ind_mfa[2]]
        lin_plot_pc_6 = vdf_pc_fit[max_ind_mfa[0], max_ind_mfa[1], :]
        lin_plot_pb_4 = vdf_pb_fit[:, max_ind_mfa[1], max_ind_mfa[2]]
        lin_plot_pb_5 = vdf_pb_fit[max_ind_mfa[0], :, max_ind_mfa[2]]
        lin_plot_pb_6 = vdf_pb_fit[max_ind_mfa[0], max_ind_mfa[1], :]
        ax = axs[1 - 1]
        img_plot = img_plot_1
        xlim = [np.min(vpar_lst), np.max(vpar_lst)]
        ylim = [np.min(vper1_lst), np.max(vper1_lst)]
        img_arr = np.ma.masked_array(img_plot, mask=np.isnan(img_plot))
        img = ax.contourf(vpar_lst, vper1_lst, img_arr, levels=levels, extend='max',
                          norm=mpl.colors.LogNorm(vmin=val_min, vmax=val_max), cmap='jet')
        ax.set_xlim(xlim[0], xlim[1])
        ax.set_ylim(ylim[0], ylim[1])
        ax.set_xlabel(r'$V_{||} (km/s)$', fontsize=15)
        ax.set_ylabel(r'$V_{\perp,1} (km/s)$', fontsize=15)
        ax.set_aspect('equal')
        ax.grid(True)

        ax = axs[2 - 1]
        img_plot = img_plot_2
        xlim = [np.min(vpar_lst), np.max(vpar_lst)]
        ylim = [np.min(vper2_lst), np.max(vper2_lst)]
        img_arr = np.ma.masked_array(img_plot, mask=np.isnan(img_plot))
        img = ax.contourf(vpar_lst, vper2_lst, img_arr, levels=levels, extend='max',
                          norm=mpl.colors.LogNorm(vmin=val_min, vmax=val_max), cmap='jet')
        ax.set_xlim(xlim[0], xlim[1])
        ax.set_ylim(ylim[0], ylim[1])
        ax.set_xlabel(r'$V_{||} (km/s)$', fontsize=15)
        ax.set_ylabel(r'$V_{\perp,2} (km/s)$', fontsize=15)
        ax.set_aspect('equal')
        ax.grid(True)

        ax = axs[3 - 1]
        img_plot = img_plot_3
        xlim = [np.min(vper1_lst), np.max(vper1_lst)]
        ylim = [np.min(vper2_lst), np.max(vper2_lst)]
        img_arr = np.ma.masked_array(img_plot, mask=np.isnan(img_plot))
        img = ax.contourf(vper1_lst, vper2_lst, img_arr, levels=levels, extend='max',
                          norm=mpl.colors.LogNorm(vmin=val_min, vmax=val_max), cmap='jet')
        ax.set_xlim(xlim[0], xlim[1])
        ax.set_ylim(ylim[0], ylim[1])
        ax.set_xlabel(r'$V_{\perp,1} (km/s)$', fontsize=15)
        ax.set_ylabel(r'$V_{\perp,2} (km/s)$', fontsize=15)
        ax.set_aspect('equal')
        ax.grid(True)
        cb_ax = fig.add_axes([.92, .35, .01, .55])
        fig.colorbar(img, orientation='vertical', cax=cb_ax, label=r'VDF ($m^{-3}\cdot s^6$)',
                     ticks=[1.e-13, 1.e-12, 1.e-11, 1.e-10, 1.e-9, 1.e-8], format=LogFormatter(10, labelOnlyBase=False))

        ax = axs[4 - 1]
        img_plot = img_plot_4
        xlim = [np.min(vpar_lst), np.max(vpar_lst)]
        ylim = [np.min(vper1_lst), np.max(vper1_lst)]
        img_arr = np.ma.masked_array(img_plot, mask=np.isnan(img_plot))
        img = ax.contourf(vpar_lst, vper1_lst, img_arr, levels=levels, extend='max',
                          norm=mpl.colors.LogNorm(vmin=val_min, vmax=val_max), cmap='jet')
        ax.set_xlim(xlim[0], xlim[1])
        ax.set_ylim(ylim[0], ylim[1])
        ax.set_xlabel(r'$V_{||} (km/s)$', fontsize=15)
        ax.set_ylabel(r'$V_{\perp,1} (km/s)$', fontsize=15)
        ax.set_aspect('equal')
        ax.grid(True)

        ax = axs[5 - 1]
        img_plot = img_plot_5
        xlim = [np.min(vpar_lst), np.max(vpar_lst)]
        ylim = [np.min(vper2_lst), np.max(vper2_lst)]
        img_arr = np.ma.masked_array(img_plot, mask=np.isnan(img_plot))
        img = ax.contourf(vpar_lst, vper2_lst, img_arr, levels=levels, extend='max',
                          norm=mpl.colors.LogNorm(vmin=val_min, vmax=val_max), cmap='jet')
        ax.set_xlim(xlim[0], xlim[1])
        ax.set_ylim(ylim[0], ylim[1])
        ax.set_xlabel(r'$V_{||} (km/s)$', fontsize=15)
        ax.set_ylabel(r'$V_{\perp,2} (km/s)$', fontsize=15)
        ax.set_aspect('equal')
        ax.grid(True)

        ax = axs[6 - 1]
        img_plot = img_plot_6
        xlim = [np.min(vper1_lst), np.max(vper1_lst)]
        ylim = [np.min(vper2_lst), np.max(vper2_lst)]
        img_arr = np.ma.masked_array(img_plot, mask=np.isnan(img_plot))
        img = ax.contourf(vper1_lst, vper2_lst, img_arr, levels=levels, extend='max',
                          norm=mpl.colors.LogNorm(vmin=val_min, vmax=val_max), cmap='jet')
        ax.set_xlim(xlim[0], xlim[1])
        ax.set_ylim(ylim[0], ylim[1])
        ax.set_xlabel(r'$V_{\perp,1} (km/s)$', fontsize=15)
        ax.set_ylabel(r'$V_{\perp,2} (km/s)$', fontsize=15)
        ax.set_aspect('equal')
        ax.grid(True)

        ax = axs[7 - 1]
        ax.scatter(vpar_lst, lin_plot_1, color='m', s=26, marker='+')
        ax.semilogy(vpar_lst, lin_plot_4, color='g')
        ax.semilogy(vpar_lst, lin_plot_pc_4, color='r')
        ax.semilogy(vpar_lst, lin_plot_pb_4, color='b')
        ax.set_xlim(np.min(vpar_lst), np.max(vpar_lst))
        ax.set_ylim(val_min, val_max)
        ax.set_xlabel(r'$V_{||} (km/s)$', fontsize=15)
        ax.set_ylabel(r'$f(v_{||})$', fontsize=15)

        ax = axs[8 - 1]
        ax.scatter(vper1_lst, lin_plot_2, color='m', s=26, marker='+')
        ax.semilogy(vper1_lst, lin_plot_5, color='g')
        ax.semilogy(vper1_lst, lin_plot_pc_5, color='r')
        ax.semilogy(vper1_lst, lin_plot_pb_5, color='b')
        ax.set_xlim(np.min(vper1_lst), np.max(vper1_lst))
        ax.set_ylim(val_min, val_max)
        ax.set_xlabel(r'$V_{\perp,1} (km/s)$', fontsize=15)
        ax.set_ylabel(r'$f(v_{\perp,1})$', fontsize=15)

        ax = axs[9 - 1]
        ax.scatter(vper2_lst, lin_plot_3, color='m', s=26, marker='+')
        ax.semilogy(vper2_lst, lin_plot_6, color='g')
        ax.semilogy(vper2_lst, lin_plot_pc_6, color='r')
        ax.semilogy(vper2_lst, lin_plot_pb_6, color='b')
        ax.set_xlim(np.min(vper2_lst), np.max(vper2_lst))
        ax.set_ylim(val_min, val_max)
        ax.set_xlabel(r'$V_{\perp,2} (km/s )$', fontsize=15)
        ax.set_ylabel(r'$f(v_{\perp,2})$', fontsize=15)
        # plt.show()
        plt.savefig(dir_fig + filename, dpi=300)
        plt.close()
