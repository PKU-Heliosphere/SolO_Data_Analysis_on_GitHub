import numpy as np
import matplotlib.pyplot as plt
from spacepy import pycdf
import os
import matplotlib as mpl
import datetime as dtime
import pycwt as wavelet
from scipy import stats
from matplotlib import gridspec
from mpl_toolkits.axes_grid1 import make_axes_locatable


def search_file(path, suffix):
    folderarray = np.empty(0)
    for file in sorted(os.listdir(path)):
        if file.endswith(suffix):
            folderarray = np.append(folderarray, file)
    return folderarray


data_type = 'pas-grnd-mom'
# acquire spherical coordinates grid
year_beg, year_end = 2021, 2021
for i_year in range(year_beg, year_end + 1):
    pas_dir = '/Users/psr/data/solar-orbiter/swa/science/l2/' + data_type + '/' + '{:04d}'.format(i_year) + '/'
    suffix = '.cdf'
    for i_file, filename in enumerate(search_file(pas_dir, suffix)):
        year, year_str = int(filename[-16:-12]), filename[-16:-12]
        mon, mon_str = int(filename[-12:-10]), filename[-12:-10]
        day, day_str = int(filename[-10:-8]), filename[-10:-8]
        # load mfi data
        pas_cdf = pycdf.CDF(pas_dir + filename)
        print(pas_cdf)
        epoch_lst = pas_cdf['Epoch'][...]
        n_lst, t_lst = pas_cdf['N'][...], pas_cdf['T'][...]
        vr_lst, vt_lst, vn_lst = pas_cdf['V_RTN'][:, 0], pas_cdf['V_RTN'][:, 1], pas_cdf['V_RTN'][:, 2]
        Tr_lst, Tt_lst, Tn_lst = pas_cdf['TxTyTz_RTN'][:, 0], pas_cdf['TxTyTz_RTN'][:, 1], pas_cdf['TxTyTz_RTN'][:, 2]
        vx_sc_lst = pas_cdf['V_SOLO_RTN'][:, 0]
        vy_sc_lst = pas_cdf['V_SOLO_RTN'][:, 1]
        vz_sc_lst = pas_cdf['V_SOLO_RTN'][:, 2]
        quality_flag = pas_cdf['quality_factor'][...]
        total_counts = pas_cdf['total_count'][...]
        helper = np.vectorize(lambda x: 1 / x.total_seconds())
        vec_time_res = helper(epoch_lst[1:] - epoch_lst[0:-1])
        vec_time_res = np.append(vec_time_res, vec_time_res[-1])
        print(stats.mode(vec_time_res))
        # quality factor: vacant!!

        # extract time series for wavelet a nalysis
        time_beg_wt, time_end_wt = dtime.datetime(2021, 11, 16, 0, 0, 0), dtime.datetime(2021, 11, 16, 12, 0, 0)
        # time_beg_wt, time_end_wt = dtime.datetime(2021, 11, 17, 0, 0, 0), dtime.datetime(2021, 11, 18, 0, 0, 0)
        wt_ind = (epoch_lst > time_beg_wt) & (epoch_lst < time_end_wt)
        epoch_wt_lst, epoch_wt_lst_v2 = epoch_lst[wt_ind], epoch_lst[wt_ind][0:-1]
        n_wt_lst = n_lst[wt_ind]
        dn_wt_lst = np.diff(n_wt_lst)
        vr_wt_lst, vt_wt_lst, vn_wt_lst = vr_lst[wt_ind], vt_lst[wt_ind], vn_lst[wt_ind]
        dvr_wt_lst, dvt_wt_lst, dvn_wt_lst = np.diff(vr_wt_lst), np.diff(vt_wt_lst), np.diff(vn_wt_lst)

        # calculate wavelet spectra using pycwt
        mwt = wavelet.Morlet(6)
        dt = float(1 / stats.mode(vec_time_res)[0])
        s0, s1 = 8, 5000  # min & max scales
        periods = np.linspace(s0, s1, 32)
        freqs = 1 / periods
        n_wave, n_scales, n_freqs, n_coi, n_fft, n_fftfreqs = wavelet.cwt(n_wt_lst, dt, wavelet=mwt, freqs=freqs)
        vr_wave, vr_scales, vr_freqs, vr_coi, vr_fft, vr_fftfreqs = wavelet.cwt(vr_wt_lst, dt, wavelet=mwt, freqs=freqs)
        vt_wave, vt_scales, vt_freqs, vt_coi, vt_fft, vt_fftfreqs = wavelet.cwt(vt_wt_lst, dt, wavelet=mwt, freqs=freqs)
        vn_wave, vn_scales, vn_freqs, vn_coi, vn_fft, vn_fftfreqs = wavelet.cwt(vn_wt_lst, dt, wavelet=mwt, freqs=freqs)
        dn_wave, dn_scales, dn_freqs, dn_coi, dn_fft, dn_fftfreqs = wavelet.cwt(dn_wt_lst, dt, wavelet=mwt, freqs=freqs)
        dvr_wave, dvr_scales, dvr_freqs, dvr_coi, dvr_fft, dvr_fftfreqs = wavelet.cwt(dvr_wt_lst, dt, wavelet=mwt,
                                                                                      freqs=freqs)
        dvt_wave, dvt_scales, dvt_freqs, dvt_coi, dvt_fft, dvt_fftfreqs = wavelet.cwt(dvt_wt_lst, dt, wavelet=mwt,
                                                                                      freqs=freqs)
        dvn_wave, dvn_scales, dvn_freqs, dvn_coi, dvn_fft, dvn_fftfreqs = wavelet.cwt(dvn_wt_lst, dt, wavelet=mwt,
                                                                                      freqs=freqs)

        # calculate psd of wavelet spectra
        n_psd_arr = np.abs(n_wave) ** 2 * 2 * dt
        vr_psd_arr = np.abs(vr_wave) ** 2 * 2 * dt
        vt_psd_arr = np.abs(vt_wave) ** 2 * 2 * dt
        vn_psd_arr = np.abs(vn_wave) ** 2 * 2 * dt
        vtrace_psd_arr = vr_psd_arr + vt_psd_arr + vn_psd_arr
        dn_psd_arr = np.abs(dn_wave) ** 2 * 2 * dt
        dvr_psd_arr = np.abs(dvr_wave) ** 2 * 2 * dt
        dvt_psd_arr = np.abs(dvt_wave) ** 2 * 2 * dt
        dvn_psd_arr = np.abs(dvn_wave) ** 2 * 2 * dt

        # calculate 1d psd of wavelet spectra
        time_beg_1dwt, time_end_1dwt = dtime.datetime(2021, 11, 16, 0, 0, 0), dtime.datetime(2021, 11, 16, 11, 0, 0)
        # time_beg_1dwt, time_end_1dwt = dtime.datetime(2021, 11, 17, 9, 0, 0), dtime.datetime(2021, 11, 17, 12, 0, 0)
        wt1d_ind = (epoch_wt_lst > time_beg_1dwt) & (epoch_wt_lst < time_end_1dwt)
        n_wt_psd_lst = np.mean(n_psd_arr[:, wt1d_ind], axis=1)
        vr_wt_psd_lst = np.mean(vr_psd_arr[:, wt1d_ind], axis=1)
        vt_wt_psd_lst = np.mean(vt_psd_arr[:, wt1d_ind], axis=1)
        vn_wt_psd_lst = np.mean(vn_psd_arr[:, wt1d_ind], axis=1)
        n_fft_psd_lst = np.abs(n_fft) ** 2 * dt * 2
        vr_fft_psd_lst = np.abs(vr_fft) ** 2 * dt * 2
        vt_fft_psd_lst = np.abs(vt_fft) ** 2 * dt * 2
        vn_fft_psd_lst = np.abs(vn_fft) ** 2 * dt * 2
        wt1d_ind = (epoch_wt_lst_v2 > time_beg_1dwt) & (epoch_wt_lst_v2 < time_end_1dwt)
        dn_wt_psd_lst = np.mean(dn_psd_arr[:, wt1d_ind], axis=1)
        dvr_wt_psd_lst = np.mean(dvr_psd_arr[:, wt1d_ind], axis=1)
        dvt_wt_psd_lst = np.mean(dvt_psd_arr[:, wt1d_ind], axis=1)
        dvn_wt_psd_lst = np.mean(dvn_psd_arr[:, wt1d_ind], axis=1)
        dn_fft_psd_lst = np.abs(dn_fft) ** 2 * dt * 2
        dvr_fft_psd_lst = np.abs(dvr_fft) ** 2 * dt * 2
        dvt_fft_psd_lst = np.abs(dvt_fft) ** 2 * dt * 2
        dvn_fft_psd_lst = np.abs(dvn_fft) ** 2 * dt * 2

        # calculate noise psd following He et al. (2019)
        # calculate structure functions
        n_sf_lst, dn_sf_lst = np.empty(len(periods)), np.empty(len(periods))
        vr_sf_lst, vt_sf_lst, vn_sf_lst = np.empty(len(periods)), np.empty(len(periods)), np.empty(len(periods))
        dvr_sf_lst, dvt_sf_lst, dvn_sf_lst = np.empty(len(periods)), np.empty(len(periods)), np.empty(len(periods))
        pixlag_lst = np.around(periods / dt / 2) * 2
        for i_v in range(8):
            if i_v == 0: var_lst = vr_wt_lst
            if i_v == 1: var_lst = vt_wt_lst
            if i_v == 2: var_lst = vn_wt_lst
            if i_v == 3: var_lst = n_wt_lst
            if i_v == 4: var_lst = dvr_wt_lst
            if i_v == 5: var_lst = dvt_wt_lst
            if i_v == 6: var_lst = dvn_wt_lst
            if i_v == 7: var_lst = dn_wt_lst
            for i_p in range(len(periods)):
                pix_shift = int(pixlag_lst[i_p])
                var_shift_lst = np.roll(var_lst, pix_shift)
                var_shift_lst[0:pix_shift - 1] = np.nan
                inc_var_lst = var_lst - var_shift_lst
                if i_v == 0: vr_sf_lst[i_p] = np.nanmean(np.abs(inc_var_lst) ** 2)
                if i_v == 1: vt_sf_lst[i_p] = np.nanmean(np.abs(inc_var_lst) ** 2)
                if i_v == 2: vn_sf_lst[i_p] = np.nanmean(np.abs(inc_var_lst) ** 2)
                if i_v == 3: n_sf_lst[i_p] = np.nanmean(np.abs(inc_var_lst) ** 2)
                if i_v == 4: dvr_sf_lst[i_p] = np.nanmean(np.abs(inc_var_lst) ** 2)
                if i_v == 5: dvt_sf_lst[i_p] = np.nanmean(np.abs(inc_var_lst) ** 2)
                if i_v == 6: dvn_sf_lst[i_p] = np.nanmean(np.abs(inc_var_lst) ** 2)
                if i_v == 6: dn_sf_lst[i_p] = np.nanmean(np.abs(inc_var_lst) ** 2)
        n_inc_sf2 = 8 / len(n_wt_lst) * n_sf_lst * 100
        vr_inc_sf2 = 8 / len(vr_wt_lst) * vr_sf_lst * 100
        vt_inc_sf2 = 8 / len(vt_wt_lst) * vt_sf_lst * 100
        vn_inc_sf2 = 8 / len(vn_wt_lst) * vn_sf_lst * 100
        n_inc_sf2 = 8 / len(n_wt_lst) * n_sf_lst * dn_wt_psd_lst / 2 / dt
        vr_inc_sf2 = 8 / len(vr_wt_lst) * vr_sf_lst * dvr_wt_psd_lst / 2 / dt
        vt_inc_sf2 = 8 / len(vt_wt_lst) * vt_sf_lst * dvt_wt_psd_lst / 2 / dt
        vn_inc_sf2 = 8 / len(vn_wt_lst) * vn_sf_lst * dvn_wt_psd_lst / 2 / dt
        n_inc_sf2 = 8 / len(n_wt_lst) * n_sf_lst * dn_sf_lst
        vr_inc_sf2 = 8 / len(vr_wt_lst) * vr_sf_lst * dvr_sf_lst
        vt_inc_sf2 = 8 / len(vt_wt_lst) * vt_sf_lst * dvt_sf_lst
        vn_inc_sf2 = 8 / len(vn_wt_lst) * vn_sf_lst * dvn_sf_lst
        alpha_n = np.mean(np.abs(dn_wt_lst - np.mean(dn_wt_lst)) ** 2) / (
                np.log10(np.max(n_freqs)) - np.log10(np.min(n_freqs)))
        alpha_vr = np.mean(np.abs(dvr_wt_lst - np.mean(dvr_wt_lst)) ** 2) / (
                np.log10(np.max(vr_freqs)) - np.log10(np.min(vr_freqs)))
        alpha_vt = np.mean(np.abs(dvt_wt_lst - np.mean(dvt_wt_lst)) ** 2) / (
                np.log10(np.max(vr_freqs)) - np.log10(np.min(vr_freqs)))
        alpha_vn = np.mean(np.abs(dvn_wt_lst - np.mean(dvn_wt_lst)) ** 2) / (
                np.log10(np.max(vr_freqs)) - np.log10(np.min(vr_freqs)))
        n_err_psd_lst = alpha_n / n_freqs * np.sqrt(n_inc_sf2)
        vr_err_psd_lst = alpha_vr / vr_freqs * np.sqrt(vr_inc_sf2)
        vt_err_psd_lst = alpha_vt / vt_freqs * np.sqrt(vt_inc_sf2)
        vn_err_psd_lst = alpha_vn / vn_freqs * np.sqrt(vn_inc_sf2)

        # plot figures
        dir_fig = '/Users/psr/work/projects/outer boundary/figures/'
        # plot time series
        file_fig = 'Time_series_solo_l2_swa-' + data_type + '_' + year_str + mon_str + day_str + '.png'
        fig, axs = plt.subplots(8, 1, figsize=[16, 16], sharex=True)
        plt.subplots_adjust(wspace=0.0, hspace=0.2, left=0.1, right=0.9, top=0.95, bottom=0.05)
        ax = axs[0]
        ax.plot(epoch_lst, vr_lst, color='k')
        ax.set_ylabel(r'$V_R$ (km/s)', fontsize=20)
        ax.tick_params(axis='y', which='major', labelsize=16)

        ax = axs[1]
        ax.plot(epoch_lst, vt_lst, color='k')
        ax.set_ylabel(r'$V_T$ (km/s)', fontsize=20)
        ax.tick_params(axis='y', which='major', labelsize=16)

        ax = axs[2]
        ax.plot(epoch_lst, vn_lst, color='k')
        ax.set_ylabel(r'$V_N$ (km/s)', fontsize=20)
        ax.tick_params(axis='y', which='major', labelsize=16)

        ax = axs[3]
        ax.plot(epoch_lst, n_lst, color='r')
        ax.set_ylabel(r'$N$ ($cm^{-3}$)', fontsize=20)
        ax.tick_params(axis='y', which='major', labelsize=16, colors='r')
        ax.spines['left'].set_color('red')
        ax2 = ax.twinx()
        ax2.plot(epoch_lst, t_lst, color='g')
        ax2.set_ylabel(r'$T$ (eV)', fontsize=20)
        ax2.tick_params(axis='y', which='major', labelsize=16, colors='g')
        ax2.spines['right'].set_color('green')

        ax = axs[4]
        ax.plot(epoch_lst, Tr_lst, color='k', label='$T_R$')
        ax.plot(epoch_lst, Tt_lst, color='r', label='$T_T$')
        ax.plot(epoch_lst, Tn_lst, color='b', label='$T_N$')
        ax.set_ylabel(r'$T$ (eV)', fontsize=20)
        ax.tick_params(axis='y', which='major', labelsize=16)
        ax.legend()

        ax = axs[5]
        ax.plot(epoch_lst, total_counts, color='k')
        ax.set_ylabel(r'Counts', fontsize=20)
        ax.tick_params(axis='y', which='major', labelsize=16)

        ax = axs[6]
        ax.plot(epoch_lst, quality_flag, color='k')
        ax.set_ylabel(r'$Q_{flag}$', fontsize=20)
        ax.tick_params(axis='y', which='major', labelsize=16)

        ax = axs[7]
        ax.plot(epoch_lst, vec_time_res, color='k')
        ax.set_ylabel(r'$f (Hz)$', fontsize=20)
        ax.tick_params(axis='y', which='major', labelsize=16)
        date_format = mpl.dates.DateFormatter('%H:%M')
        ax.xaxis.set_major_formatter(date_format)
        ax.xaxis.set_major_locator(mpl.dates.HourLocator(interval=4))
        ax.xaxis.set_minor_locator(mpl.dates.HourLocator(interval=1))
        ax.tick_params(axis='x', which='major', labelsize=16)
        fig.suptitle(year_str + '-' + mon_str + '-' + day_str, fontsize=20)
        plt.savefig(dir_fig + file_fig, dpi=300)

        # plot 2d wavelet spectra
        file_fig = 'Wavelet_Spectra_solo_l2_swa-' + data_type + '_' + year_str + mon_str + day_str + '.png'
        fig, axs = plt.subplots(4, 1, figsize=[12, 10], sharex=True)
        plt.subplots_adjust(wspace=0.0, hspace=0.25, left=0.1, right=1.00, top=0.9, bottom=0.1)
        ax = axs[0]
        vmin = np.percentile(vr_psd_arr, 5)
        vmax = np.percentile(vr_psd_arr, 95)
        levels = np.logspace(np.log10(vmin), np.log10(vmax), 20)
        im = ax.contourf(epoch_wt_lst, periods, vr_psd_arr, levels=levels, cmap='jet', extend='both',
                         norm=mpl.colors.LogNorm(vmin=vmin, vmax=vmax), shading='auto')
        ax.set_yscale('log')
        ax.set_ylabel('periods (s)', fontsize=20)
        ax.set_title('PSD($V_R$) ($(km/s)^2/Hz$)', fontsize=20)
        ax.tick_params(axis='y', which='major', labelsize=16)
        cax = plt.colorbar(im, ax=ax, format=mpl.ticker.LogFormatter(10))
        cax.ax.set_yscale('log')

        ax = axs[1]
        vmin = np.percentile(vt_psd_arr, 5)
        vmax = np.percentile(vt_psd_arr, 95)
        levels = np.logspace(np.log10(vmin), np.log10(vmax), 20)
        im = ax.contourf(epoch_wt_lst, periods, vt_psd_arr, cmap='jet', levels=levels, extend='both',
                         norm=mpl.colors.LogNorm(vmin=vmin, vmax=vmax))
        ax.set_yscale('log')
        ax.set_ylabel('periods (s)', fontsize=20)
        ax.set_title('PSD($V_T$) ($(km/s)^2/Hz$)', fontsize=20)
        ax.tick_params(axis='y', which='major', labelsize=16)
        cax = plt.colorbar(im, ax=ax, format=mpl.ticker.LogFormatter(10))
        cax.ax.set_yscale('log')

        ax = axs[2]
        vmin = np.percentile(vn_psd_arr, 5)
        vmax = np.percentile(vn_psd_arr, 95)
        levels = np.logspace(np.log10(vmin), np.log10(vmax), 20)
        im = ax.contourf(epoch_wt_lst, periods, vn_psd_arr, cmap='jet', levels=levels, extend='both',
                         norm=mpl.colors.LogNorm(vmin=vmin, vmax=vmax))
        ax.set_yscale('log')
        ax.set_ylabel('periods (s)', fontsize=20)
        ax.set_title('PSD($V_N$) ($(km/s)^2/Hz$)', fontsize=20)
        ax.tick_params(axis='y', which='major', labelsize=16)
        cax = plt.colorbar(im, ax=ax, format=mpl.ticker.LogFormatter(10))
        cax.ax.set_yscale('log')

        ax = axs[3]
        vmin = np.percentile(vtrace_psd_arr, 5)
        vmax = np.percentile(vtrace_psd_arr, 95)
        levels = np.logspace(np.log10(vmin), np.log10(vmax), 20)
        im = ax.contourf(epoch_wt_lst, periods, vtrace_psd_arr, cmap='jet', levels=levels, extend='both',
                         norm=mpl.colors.LogNorm(vmin=vmin, vmax=vmax))
        ax.set_yscale('log')
        ax.set_ylabel('periods (s)', fontsize=20)
        ax.set_title('PSD($V_{trace}$) ($(km/s)^2/Hz$)', fontsize=20)
        ax.tick_params(axis='y', which='major', labelsize=16)
        date_format = mpl.dates.DateFormatter('%H:%M')
        ax.xaxis.set_major_formatter(date_format)
        ax.xaxis.set_major_locator(mpl.dates.HourLocator(interval=3))
        ax.xaxis.set_minor_locator(mpl.dates.HourLocator(interval=1))
        ax.tick_params(axis='x', which='major', labelsize=16)
        cax = plt.colorbar(im, ax=ax, format=mpl.ticker.LogFormatter(10))
        cax.ax.set_yscale('log')
        fig.suptitle(year_str + '-' + mon_str + '-' + day_str, fontsize=20)
        plt.savefig(dir_fig + file_fig, dpi=300)

        # # plot 1d spectra
        file_fig = 'Wavelet_Spectra_1d_solo_l2_swa-' + data_type + '_' + year_str + mon_str + day_str + '.png'
        fig, axs = plt.subplots(1, 3, figsize=[16, 5], sharex=True)
        plt.subplots_adjust(wspace=0.25, hspace=0.0, left=0.1, right=0.95, top=0.9, bottom=0.1)

        ax = axs[0]
        # ax.loglog(br_fftfreqs, br_fft_psd_lst, 'g', label='fft')
        ax.loglog(vr_freqs, vr_wt_psd_lst, 'b', label='PSD')
        ax.loglog(vr_freqs, vr_err_psd_lst, 'k', label='Err for PSD')
        ax.axhline(y=1.e3, color='k')
        ax.set_ylabel(r'PSD($V_R$) ($(km/s)^2/Hz$)', fontsize=16)
        ax.set_xlabel(r'$f$ (Hz)', fontsize=16)
        ax.legend()

        ax = axs[1]
        # ax.loglog(bt_fftfreqs, bt_fft_psd_lst, 'g', label='fft')
        ax.loglog(vt_freqs, vt_wt_psd_lst, 'b', label='PSD')
        ax.loglog(vt_freqs, vt_err_psd_lst, 'k', label='Err for PSD')
        ax.axhline(y=1.e3, color='k')
        ax.set_ylabel(r'PSD($V_T$) ($(km/s)^2/Hz$)', fontsize=16)
        ax.set_xlabel(r'$f$ (Hz)', fontsize=16)
        ax.legend()

        ax = axs[2]
        # ax.loglog(bn_fftfreqs, bn_fft_psd_lst, 'g', label='fft')
        ax.loglog(vn_freqs, vn_wt_psd_lst, 'b', label='PSD')
        ax.loglog(vn_freqs, vn_err_psd_lst, 'k', label='Err for PSD')
        ax.axhline(y=1.e3, color='k')
        ax.set_ylabel(r'PSD($V_N$) ($(km/s)^2/Hz$)', fontsize=16)
        ax.set_xlabel(r'$f$ (Hz)', fontsize=16)
        ax.legend()
        plt.savefig(dir_fig + file_fig, dpi=300)

        # # plot 1d spectra
        file_fig = 'Wavelet_Spectra_dens_solo_l2_swa-' + data_type + '_' + year_str + mon_str + day_str + '.png'
        fig, ax = plt.subplots(1, 1, figsize=[10, 5])
        plt.subplots_adjust(wspace=0.2, hspace=0.2, left=0.1, right=1.00, top=0.85, bottom=0.1)
        vmin = np.percentile(n_psd_arr, 5)
        vmax = np.percentile(n_psd_arr, 95)
        levels = np.logspace(np.log10(vmin), np.log10(vmax), 20)
        im = ax.contourf(epoch_wt_lst, periods, n_psd_arr, cmap='jet', levels=levels, extend='both',
                         norm=mpl.colors.LogNorm(vmin=vmin, vmax=vmax))
        ax.set_yscale('log')
        ax.set_ylabel('periods (s)', fontsize=20)
        ax.set_title('PSD($N$) ($(cm^{-3})^2/Hz$)', fontsize=20)
        ax.tick_params(axis='y', which='major', labelsize=16)
        date_format = mpl.dates.DateFormatter('%H:%M')
        ax.xaxis.set_major_formatter(date_format)
        ax.xaxis.set_major_locator(mpl.dates.HourLocator(interval=3))
        ax.xaxis.set_minor_locator(mpl.dates.HourLocator(interval=1))
        ax.tick_params(axis='x', which='major', labelsize=16)
        cax = plt.colorbar(im, ax=ax, format=mpl.ticker.LogFormatter(10))
        cax.ax.set_yscale('log')
        fig.suptitle(year_str + '-' + mon_str + '-' + day_str, fontsize=20)
        plt.savefig(dir_fig + file_fig, dpi=300)

        # plot 1d spectra
        file_fig = 'Wavelet_Spectra_dens_1d_solo_l2_swa-' + data_type + '_' + year_str + mon_str + day_str + '.png'
        fig, ax = plt.subplots(1, 1, figsize=[6, 6], sharex=True)
        plt.subplots_adjust(wspace=0.2, hspace=0.0, left=0.15, right=0.95, top=0.9, bottom=0.1)
        ax.loglog(n_freqs, n_wt_psd_lst, 'b', label='PSD')
        ax.loglog(n_freqs, n_err_psd_lst, 'k', label='Err for PSD')
        ax.axhline(y=10, color='k')
        ax.set_ylabel(r'PSD($N$) ($(cm^{-3})^2/Hz$)', fontsize=16)
        ax.set_xlabel(r'$f$ (Hz)', fontsize=16)
        ax.legend()
        plt.savefig(dir_fig + file_fig, dpi=300)
