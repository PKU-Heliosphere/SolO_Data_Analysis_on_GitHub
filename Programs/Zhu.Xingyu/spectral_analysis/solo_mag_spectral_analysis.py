import numpy as np
import matplotlib.pyplot as plt
from spacepy import pycdf
import os
import matplotlib as mpl
import datetime as dtime
import pycwt as wavelet
from mpl_toolkits.axes_grid1 import make_axes_locatable


def search_file(path, suffix):
    folderarray = np.empty(0)
    for file in sorted(os.listdir(path)):
        if file.endswith(suffix):
            folderarray = np.append(folderarray, file)
    return folderarray


data_type = 'rtn-normal'
# acquire spherical coordinates grid
year_beg, year_end = 2021, 2021
for i_year in range(year_beg, year_end + 1):
    mag_dir = '/Users/psr/data/solar-orbiter/mag/science/l2/' + data_type + '/' + '{:04d}'.format(i_year) + '/'
    suffix = '.cdf'
    for i_file, filename in enumerate(search_file(mag_dir, suffix)):
        year, year_str = int(filename[-16:-12]), filename[-16:-12]
        mon, mon_str = int(filename[-12:-10]), filename[-12:-10]
        day, day_str = int(filename[-10:-8]), filename[-10:-8]
        # load mfi data
        mag_cdf = pycdf.CDF(mag_dir + filename)
        print(mag_cdf)
        epoch_lst = mag_cdf['EPOCH'][...]
        br_lst, bt_lst, bn_lst = mag_cdf['B_RTN'][:, 0], mag_cdf['B_RTN'][:, 1], mag_cdf['B_RTN'][:, 2]
        vec_time_res = mag_cdf['VECTOR_TIME_RESOLUTION'][...]
        quality_flag = mag_cdf['QUALITY_FLAG'][...]
        f1_lst = np.sqrt(br_lst ** 2 + bt_lst ** 2 + bn_lst ** 2)
        # quality Flag Definition
        # 0: Bad data
        # 1: Known problems use at your own risk
        # 2: Survey data, possibly not publication quality
        # 3: Good for publication, subject to PI approval
        # 4: Excellent data which has received special treatment
        # plot time series of magnetic field

        # extract time series for wavelet analysis
        time_beg_wt, time_end_wt = dtime.datetime(2021, 2, 18, 0, 0, 0), dtime.datetime(2021, 2, 18, 12, 0, 0)
        wt_ind = (epoch_lst > time_beg_wt) & (epoch_lst < time_end_wt)
        epoch_wt_lst, epoch_wt_lst_v2 = epoch_lst[wt_ind], epoch_lst[wt_ind][0:-1]
        br_wt_lst, bt_wt_lst, bn_wt_lst = br_lst[wt_ind], bt_lst[wt_ind], bn_lst[wt_ind]
        dbr_wt_lst, dbt_wt_lst, dbn_wt_lst = np.diff(br_wt_lst), np.diff(bt_wt_lst), np.diff(bn_wt_lst)

        # calculate wavelet spectra using pycwt
        mwt = wavelet.Morlet(6)
        dt = 1 / vec_time_res[0]
        s0, s1 = 0.2, 1000  # min & max scales
        periods = np.linspace(s0, s1, 32)
        freqs = 1 / periods
        br_wave, br_scales, br_freqs, br_coi, br_fft, br_fftfreqs = wavelet.cwt(br_wt_lst, dt, wavelet=mwt, freqs=freqs)
        bt_wave, bt_scales, bt_freqs, bt_coi, bt_fft, bt_fftfreqs = wavelet.cwt(bt_wt_lst, dt, wavelet=mwt, freqs=freqs)
        bn_wave, bn_scales, bn_freqs, bn_coi, bn_fft, bn_fftfreqs = wavelet.cwt(bn_wt_lst, dt, wavelet=mwt, freqs=freqs)
        dbr_wave, dbr_scales, dbr_freqs, dbr_coi, dbr_fft, dbr_fftfreqs = wavelet.cwt(dbr_wt_lst, dt, wavelet=mwt,
                                                                                      freqs=freqs)
        dbt_wave, dbt_scales, dbt_freqs, dbt_coi, dbt_fft, dbt_fftfreqs = wavelet.cwt(dbt_wt_lst, dt, wavelet=mwt,
                                                                                      freqs=freqs)
        dbn_wave, dbn_scales, dbn_freqs, dbn_coi, dbn_fft, dbn_fftfreqs = wavelet.cwt(dbn_wt_lst, dt, wavelet=mwt,
                                                                                      freqs=freqs)
        # calculate psd of wavelet spectra
        br_psd_arr = np.abs(br_wave) ** 2 * 2 * dt
        bt_psd_arr = np.abs(bt_wave) ** 2 * 2 * dt
        bn_psd_arr = np.abs(bn_wave) ** 2 * 2 * dt
        btrace_psd_arr = br_psd_arr + bt_psd_arr + bn_psd_arr
        dbr_psd_arr = np.abs(dbr_wave) ** 2 * 2 * dt
        dbt_psd_arr = np.abs(dbt_wave) ** 2 * 2 * dt
        dbn_psd_arr = np.abs(dbn_wave) ** 2 * 2 * dt

        # calculate 1d psd of wavelet spectra
        time_beg_1dwt, time_end_1dwt = dtime.datetime(2021, 2, 18, 3, 0, 0), dtime.datetime(2021, 2, 18, 9, 0, 0)
        wt1d_ind = (epoch_wt_lst > time_beg_1dwt) & (epoch_wt_lst < time_end_1dwt)
        br_wt_psd_lst = np.mean(br_psd_arr[:, wt1d_ind], axis=1)
        bt_wt_psd_lst = np.mean(bt_psd_arr[:, wt1d_ind], axis=1)
        bn_wt_psd_lst = np.mean(bn_psd_arr[:, wt1d_ind], axis=1)
        br_fft_psd_lst = np.abs(br_fft) ** 2 * dt * 2
        bt_fft_psd_lst = np.abs(bt_fft) ** 2 * dt * 2
        bn_fft_psd_lst = np.abs(bn_fft) ** 2 * dt * 2
        wt1d_ind = (epoch_wt_lst_v2 > time_beg_1dwt) & (epoch_wt_lst_v2 < time_end_1dwt)
        dbr_wt_psd_lst = np.mean(dbr_psd_arr[:, wt1d_ind], axis=1)
        dbt_wt_psd_lst = np.mean(dbt_psd_arr[:, wt1d_ind], axis=1)
        dbn_wt_psd_lst = np.mean(dbn_psd_arr[:, wt1d_ind], axis=1)
        dbr_fft_psd_lst = np.abs(dbr_fft) ** 2 * dt * 2
        dbt_fft_psd_lst = np.abs(dbt_fft) ** 2 * dt * 2
        dbn_fft_psd_lst = np.abs(dbn_fft) ** 2 * dt * 2

        # calculate noise psd following He et al. (2019)
        # calculate structure functions
        br_sf_lst, bt_sf_lst, bn_sf_lst = np.empty(len(periods)), np.empty(len(periods)), np.empty(len(periods))
        dbr_sf_lst, dbt_sf_lst, dbn_sf_lst = np.empty(len(periods)), np.empty(len(periods)), np.empty(len(periods))
        pixlag_lst = np.around(periods / dt / 2) * 2
        for i_v in range(6):
            if i_v == 0: var_lst = br_wt_lst
            if i_v == 1: var_lst = bt_wt_lst
            if i_v == 2: var_lst = bn_wt_lst
            if i_v == 3: var_lst = dbr_wt_lst
            if i_v == 4: var_lst = dbt_wt_lst
            if i_v == 5: var_lst = dbn_wt_lst
            for i_p in range(len(periods)):
                pix_shift = int(pixlag_lst[i_p])
                var_shift_lst = np.roll(var_lst, pix_shift)
                var_shift_lst[0:pix_shift - 1] = np.nan
                inc_var_lst = var_lst - var_shift_lst
                if i_v == 0: br_sf_lst[i_p] = np.nanmean(np.abs(inc_var_lst) ** 2)
                if i_v == 1: bt_sf_lst[i_p] = np.nanmean(np.abs(inc_var_lst) ** 2)
                if i_v == 2: bn_sf_lst[i_p] = np.nanmean(np.abs(inc_var_lst) ** 2)
                if i_v == 3: dbr_sf_lst[i_p] = np.nanmean(np.abs(inc_var_lst) ** 2)
                if i_v == 4: dbt_sf_lst[i_p] = np.nanmean(np.abs(inc_var_lst) ** 2)
                if i_v == 5: dbn_sf_lst[i_p] = np.nanmean(np.abs(inc_var_lst) ** 2)
        br_inc_sf2 = 8 / len(br_wt_lst) * br_sf_lst * 0.01
        bt_inc_sf2 = 8 / len(bt_wt_lst) * bt_sf_lst * 0.01
        bn_inc_sf2 = 8 / len(bn_wt_lst) * bn_sf_lst * 0.01
        br_inc_sf2 = 8 / len(br_wt_lst) * br_sf_lst * dbr_wt_psd_lst / 2 / dt
        bt_inc_sf2 = 8 / len(bt_wt_lst) * bt_sf_lst * dbt_wt_psd_lst / 2 / dt
        bn_inc_sf2 = 8 / len(bn_wt_lst) * bn_sf_lst * dbn_wt_psd_lst / 2 / dt
        br_inc_sf2 = 8 / len(br_wt_lst) * br_sf_lst * dbr_sf_lst
        bt_inc_sf2 = 8 / len(bt_wt_lst) * bt_sf_lst * dbt_sf_lst
        bn_inc_sf2 = 8 / len(bn_wt_lst) * bn_sf_lst * dbn_sf_lst
        alpha_br = np.mean(np.abs(dbr_wt_lst - np.mean(dbr_wt_lst)) ** 2) / (
                np.log10(np.max(br_freqs)) - np.log10(np.min(br_freqs)))
        alpha_bt = np.mean(np.abs(dbt_wt_lst - np.mean(dbt_wt_lst)) ** 2) / (
                np.log10(np.max(br_freqs)) - np.log10(np.min(br_freqs)))
        alpha_bn = np.mean(np.abs(dbn_wt_lst - np.mean(dbn_wt_lst)) ** 2) / (
                np.log10(np.max(br_freqs)) - np.log10(np.min(br_freqs)))
        br_err_psd_lst = alpha_br / br_freqs * np.sqrt(br_inc_sf2)
        bt_err_psd_lst = alpha_bt / bt_freqs * np.sqrt(bt_inc_sf2)
        bn_err_psd_lst = alpha_bn / bn_freqs * np.sqrt(bn_inc_sf2)

        # plot figures
        dir_fig = '/Users/psr/work/projects/outer boundary/figures/'
        # plot time series
        # file_fig = 'Time_series_solo_l2_mag-' + data_type + '_' + year_str + mon_str + day_str + '.png'
        # fig, axs = plt.subplots(6, 1, figsize=[12, 10], sharex=True)
        # plt.subplots_adjust(wspace=0.0, hspace=0.2, left=0.1, right=0.9, top=0.95, bottom=0.1)
        # ax = axs[0]
        # ax.plot(epoch_lst, f1_lst, color='k')
        # ax.set_ylabel('|B| (nT)', fontsize=20)
        # ax.tick_params(axis='y', which='major', labelsize=16)
        #
        # ax = axs[1]
        # ax.plot(epoch_lst, br_lst, color='k')
        # ax.set_ylabel(r'$B_R$ (nT)', fontsize=20)
        # ax.tick_params(axis='y', which='major', labelsize=16)
        #
        # ax = axs[2]
        # ax.plot(epoch_lst, bt_lst, color='k')
        # ax.set_ylabel(r'$B_T$ (nT)', fontsize=20)
        # ax.tick_params(axis='y', which='major', labelsize=16)
        #
        # ax = axs[3]
        # ax.plot(epoch_lst, bn_lst, color='k')
        # ax.set_ylabel(r'$B_N$ (nT)', fontsize=20)
        #
        # ax = axs[4]
        # ax.plot(epoch_lst, quality_flag, color='k')
        # ax.set_ylabel(r'$Q_{flag}$', fontsize=20)
        # ax.tick_params(axis='y', which='major', labelsize=16)
        #
        # ax = axs[5]
        # ax.plot(epoch_lst, vec_time_res, color='k')
        # ax.set_ylabel(r'$f$ (Hz)', fontsize=20)
        # ax.tick_params(axis='y', which='major', labelsize=16)
        # date_format = mpl.dates.DateFormatter('%H:%M')
        # ax.xaxis.set_major_formatter(date_format)
        # ax.xaxis.set_major_locator(mpl.dates.HourLocator(interval=4))
        # ax.xaxis.set_minor_locator(mpl.dates.HourLocator(interval=1))
        # ax.tick_params(axis='x', which='major', labelsize=16)
        # fig.suptitle(year_str + '-' + mon_str + '-' + day_str, fontsize=20)
        # plt.savefig(dir_fig + file_fig, dpi=300)
        #
        # # plot 2d wavelet spectra
        # file_fig = 'Wavelet_Spectra_solo_l2_mag-' + data_type + '_' + year_str + mon_str + day_str + '.png'
        # fig, axs = plt.subplots(4, 1, figsize=[12, 10], sharex=True)
        # plt.subplots_adjust(wspace=0.0, hspace=0.25, left=0.1, right=1.00, top=0.9, bottom=0.1)
        # ax = axs[0]
        # vmin = np.percentile(br_psd_arr, 5)
        # vmax = np.percentile(br_psd_arr, 95)
        # levels = np.logspace(np.log10(vmin), np.log10(vmax), 20)
        # im = ax.contourf(epoch_wt_lst, periods, br_psd_arr, levels=levels, cmap='jet', extend='both',
        #                  norm=mpl.colors.LogNorm(vmin=vmin, vmax=vmax), shading='auto')
        # ax.set_yscale('log')
        # ax.set_ylabel('periods (s)', fontsize=20)
        # ax.set_title('PSD($B_R$) ($nT^2/Hz$)', fontsize=20)
        # ax.tick_params(axis='y', which='major', labelsize=16)
        # cax = plt.colorbar(im, ax=ax, format=mpl.ticker.LogFormatter(10))
        # cax.ax.set_yscale('log')
        #
        # ax = axs[1]
        # vmin = np.percentile(bt_psd_arr, 5)
        # vmax = np.percentile(bt_psd_arr, 95)
        # levels = np.logspace(np.log10(vmin), np.log10(vmax), 20)
        # im = ax.contourf(epoch_wt_lst, periods, bt_psd_arr, cmap='jet', levels=levels, extend='both',
        #                  norm=mpl.colors.LogNorm(vmin=vmin, vmax=vmax))
        # ax.set_yscale('log')
        # ax.set_ylabel('periods (s)', fontsize=20)
        # ax.set_title('PSD($B_T$) ($nT^2/Hz$)', fontsize=20)
        # ax.tick_params(axis='y', which='major', labelsize=16)
        # cax = plt.colorbar(im, ax=ax, format=mpl.ticker.LogFormatter(10))
        # cax.ax.set_yscale('log')
        #
        # ax = axs[2]
        # vmin = np.percentile(bn_psd_arr, 5)
        # vmax = np.percentile(bn_psd_arr, 95)
        # levels = np.logspace(np.log10(vmin), np.log10(vmax), 20)
        # im = ax.contourf(epoch_wt_lst, periods, bn_psd_arr, cmap='jet', levels=levels, extend='both',
        #                  norm=mpl.colors.LogNorm(vmin=vmin, vmax=vmax))
        # ax.set_yscale('log')
        # ax.set_ylabel('periods (s)', fontsize=20)
        # ax.set_title('PSD($B_N$) ($nT^2/Hz$)', fontsize=20)
        # ax.tick_params(axis='y', which='major', labelsize=16)
        # cax = plt.colorbar(im, ax=ax, format=mpl.ticker.LogFormatter(10))
        # cax.ax.set_yscale('log')
        #
        # ax = axs[3]
        # vmin = np.percentile(btrace_psd_arr, 5)
        # vmax = np.percentile(btrace_psd_arr, 95)
        # levels = np.logspace(np.log10(vmin), np.log10(vmax), 20)
        # im = ax.contourf(epoch_wt_lst, periods, btrace_psd_arr, cmap='jet', levels=levels, extend='both',
        #                  norm=mpl.colors.LogNorm(vmin=vmin, vmax=vmax))
        # ax.set_yscale('log')
        # ax.set_ylabel('periods (s)', fontsize=20)
        # ax.set_title('PSD($B_{trace}$) ( $nT^2/Hz$)', fontsize=20)
        # ax.tick_params(axis='y', which='major', labelsize=16)
        # date_format = mpl.dates.DateFormatter('%H:%M')
        # ax.xaxis.set_major_formatter(date_format)
        # ax.xaxis.set_major_locator(mpl.dates.HourLocator(interval=3))
        # ax.xaxis.set_minor_locator(mpl.dates.HourLocator(interval=1))
        # ax.tick_params(axis='x', which='major', labelsize=16)
        # cax = plt.colorbar(im, ax=ax, format=mpl.ticker.LogFormatter(10))
        # cax.ax.set_yscale('log')
        # fig.suptitle(year_str + '-' + mon_str + '-' + day_str, fontsize=20)
        # plt.savefig(dir_fig + file_fig, dpi=300)

        # plot 1d spectra
        file_fig = 'Wavelet_Spectra_1d_solo_l2_mag-' + data_type + '_' + year_str + mon_str + day_str + '.png'
        fig, axs = plt.subplots(1, 3, figsize=[16, 5], sharex=True)
        plt.subplots_adjust(wspace=0.25, hspace=0.0, left=0.1, right=0.95, top=0.9, bottom=0.1)

        ax = axs[0]
        # ax.loglog(br_fftfreqs, br_fft_psd_lst, 'g', label='fft')
        ax.loglog(br_freqs, br_wt_psd_lst, 'b', label='PSD')
        ax.loglog(br_freqs, br_err_psd_lst, 'k', label='Err for PSD')
        ax.axhline(y=1.e-4, color='k')
        ax.set_ylabel(r'PSD($B_R$) ($nT^2/Hz$)', fontsize=16)
        ax.set_xlabel(r'$f$ (Hz)', fontsize=16)
        ax.legend()

        ax = axs[1]
        # ax.loglog(bt_fftfreqs, bt_fft_psd_lst, 'g', label='fft')
        ax.loglog(bt_freqs, bt_wt_psd_lst, 'b', label='PSD')
        ax.loglog(bt_freqs, bt_err_psd_lst, 'k', label='Err for PSD')
        ax.axhline(y=1.e-4, color='k')
        ax.set_ylabel(r'PSD($B_T$) ($nT^2/Hz$)', fontsize=16)
        ax.set_xlabel(r'$f$ (Hz)', fontsize=16)
        ax.legend()

        ax = axs[2]
        # ax.loglog(bn_fftfreqs, bn_fft_psd_lst, 'g', label='fft')
        ax.loglog(bn_freqs, bn_wt_psd_lst, 'b', label='PSD')
        ax.loglog(bn_freqs, bn_err_psd_lst, 'k', label='Err for PSD')
        ax.axhline(y=1.e-4, color='k')
        ax.set_ylabel(r'PSD($B_N$) ($nT^2/Hz$)', fontsize=16)
        ax.set_xlabel(r'$f$ (Hz)', fontsize=16)
        ax.legend()
        plt.savefig(dir_fig + file_fig, dpi=300)
