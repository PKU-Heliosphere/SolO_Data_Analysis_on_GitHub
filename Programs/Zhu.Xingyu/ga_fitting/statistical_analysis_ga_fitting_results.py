import numpy as np
import os
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import datetime as dt
from scipy.optimize import curve_fit
import bisect
import scipy.constants as const
import datetime as dt
from spacepy import pycdf
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter, AutoMinorLocator)


def search_file(path, suffix):
    folderarray = np.empty(0)
    for file in sorted(os.listdir(path)):
        if file.endswith(suffix):
            folderarray = np.append(folderarray, file)
    return folderarray


def plot_confidence_level():
    print('========99% confidence level========')
    print('σ_dens_pc_low, σ_dens_pc_up:', np.mean(delta_dens_pc_low_lst), ',', np.mean(delta_dens_pc_up_lst))
    print('σ_dens_pb_low, σ_dens_pb_up:', np.mean(delta_dens_pb_low_lst), ',', np.mean(delta_dens_pb_up_lst))
    print('σ_wthp_pc_low, σ_wthp_pc_up:', np.mean(delta_wthp_pc_low_lst), ',', np.mean(delta_wthp_pc_up_lst))
    print('σ_wthp_pb_low, σ_wthp_pb_up:', np.mean(delta_wthp_pb_low_lst), ',', np.mean(delta_wthp_pb_up_lst))
    print('σ_wthz_pc_low, σ_wthz_pc_up:', np.mean(delta_wthz_pc_low_lst), ',', np.mean(delta_wthz_pc_up_lst))
    print('σ_wthz_pb_low, σ_wthz_pb_up:', np.mean(delta_wthz_pb_low_lst), ',', np.mean(delta_wthz_pb_up_lst))
    print('σ_ux_pc_low, σ_ux_pc_up:', np.mean(delta_ux_pc_low_lst), ',', np.mean(delta_ux_pc_up_lst))
    print('σ_ux_pb_low, σ_ux_pb_up:', np.mean(delta_ux_pb_low_lst), ',', np.mean(delta_ux_pb_up_lst))
    print('σ_uy_pc_low, σ_uy_pc_up:', np.mean(delta_uy_pc_low_lst), ',', np.mean(delta_uy_pc_up_lst))
    print('σ_uy_pb_low, σ_uy_pb_up:', np.mean(delta_uy_pb_low_lst), ',', np.mean(delta_uy_pb_up_lst))
    print('σ_uz_pc_low, σ_uz_pc_up:', np.mean(delta_uz_pc_low_lst), ',', np.mean(delta_uz_pc_up_lst))
    print('σ_uz_pb_low, σ_uz_pb_up:', np.mean(delta_uz_pb_low_lst), ',', np.mean(delta_uz_pb_up_lst))
    print('σ_ur_pc_low, σ_ur_pc_up:', np.mean(delta_ur_pc_low_lst), ',', np.mean(delta_ur_pc_up_lst))
    print('σ_ur_pb_low, σ_ur_pb_up:', np.mean(delta_ur_pb_low_lst), ',', np.mean(delta_ur_pb_up_lst))
    print('σ_ut_pc_low, σ_ut_pc_up:', np.mean(delta_ut_pc_low_lst), ',', np.mean(delta_ut_pc_up_lst))
    print('σ_ut_pb_low, σ_ut_pb_up:', np.mean(delta_ut_pb_low_lst), ',', np.mean(delta_ut_pb_up_lst))
    print('σ_un_pc_low, σ_un_pc_up:', np.mean(delta_un_pc_low_lst), ',', np.mean(delta_un_pc_up_lst))
    print('σ_un_pb_low, σ_un_pb_up:', np.mean(delta_un_pb_low_lst), ',', np.mean(delta_un_pb_up_lst))


# set calculation time
year, month, day = 2020, 8, 2
_year, _month, _day = '{:04d}'.format(year), '{:02d}'.format(month), '{:02d}'.format(day)
beg_time = dt.datetime(year, month, day, 11, 30, 00, 000000)
end_time = dt.datetime(year, month, day, 12, 00, 00, 000000)
_date = beg_time.strftime('%Y%m%d')
_beg_t, _end_t = beg_time.strftime('%H%M%S'), end_time.strftime('%H%M%S')

# create parameter lists
datetime_vect = []
hour_vect, min_vect, sec_vect = [], [], []
# fitting parameters
dens_pc_vect, dens_pb_vect = [], []
wthp_pc_vect, wthz_pc_vect, wthp_pb_vect, wthz_pb_vect = [], [], [], []
ux_pc_vect, uy_pc_vect, uz_pc_vect = [], [], []
ux_pb_vect, uy_pb_vect, uz_pb_vect = [], [], []
file_dir = '/Users/psr/work/pku/nfa_beam/Data/solutions/' + '-'.join(
    ['{:04d}'.format(year), '{:02d}'.format(month), '{:02d}'.format(day)]) + '/'
for i_file, filename in enumerate(search_file(file_dir, 'txt')):
    _hour, _minute, _second = filename[-14:-12], filename[-12:-10], filename[-10:-4]
    dt_tmp = dt.datetime.strptime(' '.join([_year, _month, _day, _hour, _minute, _second]), '%Y %m %d %H %M %S.%f')
    dens_pc, dens_pb, wthp_pc, wthz_pc, wthp_pb, wthz_pb, ux_pc, uy_pc, uz_pc, ux_pb, uy_pb, uz_pb = np.loadtxt(
        file_dir + filename)
    dens_pc_vect.append(dens_pc)
    dens_pb_vect.append(dens_pb)
    wthp_pc_vect.append(wthp_pc)
    wthz_pc_vect.append(wthz_pc)
    wthp_pb_vect.append(wthp_pb)
    wthz_pb_vect.append(wthz_pb)
    ux_pc_vect.append(ux_pc)
    ux_pb_vect.append(ux_pb)
    uy_pc_vect.append(uy_pc)
    uy_pb_vect.append(uy_pb)
    uz_pc_vect.append(uz_pc)
    uz_pb_vect.append(uz_pb)
    datetime_vect.append(dt_tmp)
    hour_vect.append(int(_hour))
    min_vect.append(int(_minute))
    sec_vect.append(float(_second))
beg_ind = bisect.bisect_left(datetime_vect, beg_time)
end_ind = bisect.bisect_right(datetime_vect, end_time)
datetime_lst = datetime_vect[beg_ind:end_ind]
hour_lst, min_lst, sec_lst = hour_vect[beg_ind:end_ind], min_vect[beg_ind:end_ind], sec_vect[beg_ind:end_ind]
dens_pc_lst, dens_pb_lst = dens_pc_vect[beg_ind:end_ind], dens_pb_vect[beg_ind:end_ind]
wthp_pc_lst, wthp_pb_lst = wthp_pc_vect[beg_ind:end_ind], wthp_pb_vect[beg_ind:end_ind]
wthz_pc_lst, wthz_pb_lst = wthz_pc_vect[beg_ind:end_ind], wthz_pb_vect[beg_ind:end_ind]
ux_pc_lst, ux_pb_lst = ux_pc_vect[beg_ind:end_ind], ux_pb_vect[beg_ind:end_ind]
uy_pc_lst, uy_pb_lst = uy_pc_vect[beg_ind:end_ind], uy_pb_vect[beg_ind:end_ind]
uz_pc_lst, uz_pb_lst = uz_pc_vect[beg_ind:end_ind], uz_pb_vect[beg_ind:end_ind]

# read proton moments cdf file
dir_cdf = '/Users/psr/work/pku/nfa_beam/Data/'
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

# fitted bulk velocities plus maximum velocity
ur_pc_lst, ut_pc_lst, un_pc_lst = [], [], []
ur_pb_lst, ut_pb_lst, un_pb_lst = [], [], []
# lower and upper limits of fitting parameters (mfa coordinates)
dens_pc_low_lst, dens_pb_low_lst, dens_pc_up_lst, dens_pb_up_lst = [], [], [], []
wthp_pc_low_lst, wthp_pb_low_lst, wthp_pc_up_lst, wthp_pb_up_lst = [], [], [], []
wthz_pc_low_lst, wthz_pb_low_lst, wthz_pc_up_lst, wthz_pb_up_lst = [], [], [], []
ux_pc_low_lst, ux_pb_low_lst, ux_pc_up_lst, ux_pb_up_lst = [], [], [], []
uy_pc_low_lst, uy_pb_low_lst, uy_pc_up_lst, uy_pb_up_lst = [], [], [], []
uz_pc_low_lst, uz_pb_low_lst, uz_pc_up_lst, uz_pb_up_lst = [], [], [], []
# lower and upper limits of fitting parameters (bulk velocities in rtn coordinates)
ur_pc_low_lst, ur_pb_low_lst, ur_pc_up_lst, ur_pb_up_lst = [], [], [], []
ut_pc_low_lst, ut_pb_low_lst, ut_pc_up_lst, ut_pb_up_lst = [], [], [], []
un_pc_low_lst, un_pb_low_lst, un_pc_up_lst, un_pb_up_lst = [], [], [], []
# lower and upper parts of deviations of fitted parameters
delta_dens_pc_low_lst, delta_dens_pb_low_lst, delta_dens_pc_up_lst, delta_dens_pb_up_lst = [], [], [], []
delta_wthp_pc_low_lst, delta_wthp_pb_low_lst, delta_wthp_pc_up_lst, delta_wthp_pb_up_lst = [], [], [], []
delta_wthz_pc_low_lst, delta_wthz_pb_low_lst, delta_wthz_pc_up_lst, delta_wthz_pb_up_lst = [], [], [], []
delta_ux_pc_low_lst, delta_ux_pb_low_lst, delta_ux_pc_up_lst, delta_ux_pb_up_lst = [], [], [], []
delta_uy_pc_low_lst, delta_uy_pb_low_lst, delta_uy_pc_up_lst, delta_uy_pb_up_lst = [], [], [], []
delta_uz_pc_low_lst, delta_uz_pb_low_lst, delta_uz_pc_up_lst, delta_uz_pb_up_lst = [], [], [], []
delta_ur_pc_low_lst, delta_ur_pb_low_lst, delta_ur_pc_up_lst, delta_ur_pb_up_lst = [], [], [], []
delta_ut_pc_low_lst, delta_ut_pb_low_lst, delta_ut_pc_up_lst, delta_ut_pb_up_lst = [], [], [], []
delta_un_pc_low_lst, delta_un_pb_low_lst, delta_un_pc_up_lst, delta_un_pb_up_lst = [], [], [], []
# theta arrays
theta_mean_lst, theta_std_lst, theta_fit_lst, theta_ste_lst = [], [], [], []
theta_std_lst_low, theta_std_lst_up, theta_ste_lst_low, theta_ste_lst_up = [], [], [], []
# whether use Δχ2 around minimum χ2
is_find_chi_square_min = 1
for i_time, epoch in enumerate(epoch_lst):
    vdf = vdf_lst[i_time, :, :, :]
    if dt.datetime(2020, 8, 2, 11, 38, 34) < epoch < dt.datetime(2020, 8, 2, 11, 38, 35):
        vdf[np.where(vel_arr < 310)] = 0.0
    else:
        vdf[np.where(vel_arr < 300)] = 0.0
    vdf[np.where(vdf == 0.0)] = np.nan
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
    rtn2mfa = np.vstack([epara, eperp1, eperp2])
    mfa2rtn = np.linalg.inv(rtn2mfa)
    # find the maximum vdf in the original array
    if epoch < dt.datetime(2020, 8, 2, 11, 32, 0) or (
            dt.datetime(2020, 8, 2, 11, 42, 10) < epoch < dt.datetime(2020, 8, 2, 11, 42, 20)) or (
            dt.datetime(2020, 8, 2, 11, 47, 14) < epoch < dt.datetime(2020, 8, 2, 11, 47, 15)) or (
            dt.datetime(2020, 8, 2, 11, 54, 14) < epoch < dt.datetime(2020, 8, 2, 11, 54, 15)) or (
            dt.datetime(2020, 8, 2, 11, 57, 58) < epoch < dt.datetime(2020, 8, 2, 11, 57, 59)) or (
            dt.datetime(2020, 8, 2, 11, 58, 18) < epoch < dt.datetime(2020, 8, 2, 11, 58, 19)):
        vdf[np.where(
            np.logical_and(np.logical_or(np.abs(vperp1_rtn) > 150, np.abs(vperp2_rtn > 150)), vdf > 1.e-9))] = np.nan
    max_ind = np.nanargmax(vdf)
    max_vdf_tmp = np.nanmax(vdf[np.logical_and(vel_arr > 300, vel_arr < 350)])
    if dt.datetime(2020, 8, 2, 11, 38, 34) < epoch < dt.datetime(2020, 8, 2, 11, 38, 35):
        max_vdf_tmp = np.nanmax(vdf[np.logical_and(vel_arr > 310, vel_arr < 350)])
    max_ind = np.ravel_multi_index(np.where(vdf == max_vdf_tmp), (11, 9, 96))
    if dt.datetime(2020, 8, 2, 11, 57, 58) < epoch < dt.datetime(2020, 8, 2, 11, 57, 59): max_ind = 3997
    vpara_max = vpara_rtn[np.unravel_index(max_ind, (11, 9, 96))]
    vperp1_max = vperp1_rtn[np.unravel_index(max_ind, (11, 9, 96))]
    vperp2_max = vperp2_rtn[np.unravel_index(max_ind, (11, 9, 96))]
    vpara = vpara_rtn - vpara_max
    vperp1 = vperp1_rtn - vperp1_max
    vperp2 = vperp2_rtn - vperp2_max
    # calculate fitted bulk velocities in rtn coordinates
    u_pc_vect = np.vstack(
        [ux_pc_lst[i_time] + vpara_max, uy_pc_lst[i_time] + vperp1_max, uz_pc_lst[i_time] + vperp2_max])
    u_pb_vect = np.vstack(
        [ux_pb_lst[i_time] + vpara_max, uy_pb_lst[i_time] + vperp1_max, uz_pb_lst[i_time] + vperp2_max])
    u_rtn_pc_vect = np.dot(mfa2rtn, u_pc_vect)
    u_rtn_pb_vect = np.dot(mfa2rtn, u_pb_vect)
    ur_pc_lst.append(u_rtn_pc_vect[0])
    ut_pc_lst.append(u_rtn_pc_vect[1])
    un_pc_lst.append(u_rtn_pc_vect[2])
    ur_pb_lst.append(u_rtn_pb_vect[0])
    ut_pb_lst.append(u_rtn_pb_vect[1])
    un_pb_lst.append(u_rtn_pb_vect[2])
    # calculate X2 for individual parameters
    num_var = 200
    var_n_pc_lst, var_n_pb_lst = np.linspace(0, 10, num_var), np.linspace(0, 10, num_var)
    var_wthp_pc_lst, var_wthp_pb_lst = np.linspace(20, 80, num_var), np.linspace(40, 100, num_var)
    var_wthz_pc_lst, var_wthz_pb_lst = np.linspace(15, 60, num_var), np.linspace(40, 100, num_var)
    var_ux_pc_lst, var_ux_pb_lst = np.linspace(-20, 20, num_var), np.linspace(-100, -20, num_var)
    var_uy_pc_lst, var_uy_pb_lst = np.linspace(-50, 50, num_var), np.linspace(-50, 50, num_var)
    var_uz_pc_lst, var_uz_pb_lst = np.linspace(-50, 50, num_var), np.linspace(-50, 50, num_var)
    chi_n_pc_lst, chi_n_pb_lst = [], []
    chi_wthp_pc_lst, chi_wthp_pb_lst = [], []
    chi_wthz_pc_lst, chi_wthz_pb_lst = [], []
    chi_ux_pc_lst, chi_ux_pb_lst = [], []
    chi_uy_pc_lst, chi_uy_pb_lst = [], []
    chi_uz_pc_lst, chi_uz_pb_lst = [], []
    # calculate X2 for individual parameters (rtn coordinates)
    var_ur_pc_lst, var_ur_pb_lst = np.linspace(ur_pc_lst[i_time] - 35, ur_pc_lst[i_time] + 35, num_var), np.linspace(
        ur_pb_lst[i_time] - 35, ur_pb_lst[i_time] + 35, num_var)
    var_ut_pc_lst, var_ut_pb_lst = np.linspace(ut_pc_lst[i_time] - 35, ut_pc_lst[i_time] + 35, num_var), np.linspace(
        ut_pb_lst[i_time] - 35, ut_pb_lst[i_time] + 35, num_var)
    var_un_pc_lst, var_un_pb_lst = np.linspace(un_pc_lst[i_time] - 35, un_pc_lst[i_time] + 35, num_var), np.linspace(
        un_pb_lst[i_time] - 35, un_pb_lst[i_time] + 35, num_var)
    chi_ur_pc_lst, chi_ur_pb_lst = [], []
    chi_ut_pc_lst, chi_ut_pb_lst = [], []
    chi_un_pc_lst, chi_un_pb_lst = [], []
    # calculate original vdf_fit with fitted parameters
    vdf_fit = 1.e-3 * dens_pc_lst[i_time] / (wthp_pc_lst[i_time] ** 2 * wthz_pc_lst[i_time] * np.pi ** 1.5) * \
              np.exp(-(vpara - ux_pc_lst[i_time]) ** 2 / wthz_pc_lst[i_time] ** 2) * \
              np.exp(-(vperp1 - uy_pc_lst[i_time]) ** 2 / wthp_pc_lst[i_time] ** 2) * \
              np.exp(-(vperp2 - uz_pc_lst[i_time]) ** 2 / wthp_pc_lst[i_time] ** 2) + \
              1.e-3 * dens_pb_lst[i_time] / (wthp_pb_lst[i_time] ** 2 * wthz_pb_lst[i_time] * np.pi ** 1.5) * \
              np.exp(-(vpara - ux_pb_lst[i_time]) ** 2 / wthz_pb_lst[i_time] ** 2) * \
              np.exp(-(vperp1 - uy_pb_lst[i_time]) ** 2 / wthp_pb_lst[i_time] ** 2) * \
              np.exp(-(vperp2 - uz_pb_lst[i_time]) ** 2 / wthp_pb_lst[i_time] ** 2)
    sigma2_aver = np.nansum((vdf - vdf_fit) ** 2 / (np.count_nonzero(~np.isnan(vdf)) - 12))
    sigma = sigma2_aver
    mean_times_ori_sd = 0

    # def Gauss(x, a, x0, sigma):
    #     return a * np.exp(-(x - x0) ** 2 / (2 * sigma ** 2))
    #
    #
    # bins = np.linspace(-0.5, 0.5, 30)
    # bin_centers = 0.5 * (bins[1:] + bins[:-1])
    # histogram, bins = np.histogram((vdf.reshape(-1) - vdf_fit.reshape(-1)) / np.sqrt(sigma2_aver), bins=bins,
    #                                density=True)
    # mean = sum(bin_centers * histogram) / sum(histogram)
    # sigma = np.sqrt(sum(histogram * (bin_centers - mean) ** 2) / sum(histogram))
    # popt, pcov = curve_fit(Gauss, bin_centers, histogram, p0=[max(histogram), mean, sigma])

    # calculate effective sigma
    # sigma2 = popt[2] ** 2 * sigma2_aver
    # bins = np.linspace(-3, 3, 30)
    # bin_centers = 0.5 * (bins[1:] + bins[:-1])
    # histogram, bins = np.histogram(
    #     (vdf.reshape(-1) - vdf_fit.reshape(-1) - popt[1] * np.sqrt(sigma2_aver)) / np.sqrt(sigma2), bins=bins,
    #     density=True)
    # mean_times_ori_sd = popt[1] * np.sqrt(sigma2_aver)

    # calculate Δχ2 of individual parameter with other parameters fixed (bulk velocities in rtn coordinates)
    for i_chi in range(num_var):
        ux_tmp, uy_tmp, uz_tmp = np.dot(rtn2mfa, np.vstack(
            [var_ur_pc_lst[i_chi], ut_pc_lst[i_time], un_pc_lst[i_time]])) - np.vstack(
            [vpara_max, vperp1_max, vperp2_max])
        vdf_fit_ur_pc_tmp = 1.e-3 * dens_pc_lst[i_time] / (
                wthp_pc_lst[i_time] ** 2 * wthz_pc_lst[i_time] * np.pi ** 1.5) * \
                            np.exp(-(vpara - ux_tmp) ** 2 / wthz_pc_lst[i_time] ** 2) * \
                            np.exp(-(vperp1 - uy_tmp) ** 2 / wthp_pc_lst[i_time] ** 2) * \
                            np.exp(-(vperp2 - uz_tmp) ** 2 / wthp_pc_lst[i_time] ** 2) + \
                            1.e-3 * dens_pb_lst[i_time] / (
                                    wthp_pb_lst[i_time] ** 2 * wthz_pb_lst[i_time] * np.pi ** 1.5) * \
                            np.exp(-(vpara - ux_pb_lst[i_time]) ** 2 / wthz_pb_lst[i_time] ** 2) * \
                            np.exp(-(vperp1 - uy_pb_lst[i_time]) ** 2 / wthp_pb_lst[i_time] ** 2) * \
                            np.exp(-(vperp2 - uz_pb_lst[i_time]) ** 2 / wthp_pb_lst[i_time] ** 2)
        chi_ur_pc_lst.append(np.nansum((vdf - vdf_fit_ur_pc_tmp) ** 2 / sigma))
    for i_chi in range(num_var):
        ux_tmp, uy_tmp, uz_tmp = np.dot(rtn2mfa, np.vstack(
            [ur_pc_lst[i_time], var_ut_pc_lst[i_chi], un_pc_lst[i_time]])) - np.vstack(
            [vpara_max, vperp1_max, vperp2_max])
        vdf_fit_ut_pc_tmp = 1.e-3 * dens_pc_lst[i_time] / (
                wthp_pc_lst[i_time] ** 2 * wthz_pc_lst[i_time] * np.pi ** 1.5) * \
                            np.exp(-(vpara - ux_tmp) ** 2 / wthz_pc_lst[i_time] ** 2) * \
                            np.exp(-(vperp1 - uy_tmp) ** 2 / wthp_pc_lst[i_time] ** 2) * \
                            np.exp(-(vperp2 - uz_tmp) ** 2 / wthp_pc_lst[i_time] ** 2) + \
                            1.e-3 * dens_pb_lst[i_time] / (
                                    wthp_pb_lst[i_time] ** 2 * wthz_pb_lst[i_time] * np.pi ** 1.5) * \
                            np.exp(-(vpara - ux_pb_lst[i_time]) ** 2 / wthz_pb_lst[i_time] ** 2) * \
                            np.exp(-(vperp1 - uy_pb_lst[i_time]) ** 2 / wthp_pb_lst[i_time] ** 2) * \
                            np.exp(-(vperp2 - uz_pb_lst[i_time]) ** 2 / wthp_pb_lst[i_time] ** 2)
        chi_ut_pc_lst.append(np.nansum((vdf - vdf_fit_ut_pc_tmp) ** 2 / sigma))
    for i_chi in range(num_var):
        ux_tmp, uy_tmp, uz_tmp = np.dot(rtn2mfa, np.vstack(
            [ur_pc_lst[i_time], ut_pc_lst[i_time], var_un_pc_lst[i_chi]])) - np.vstack(
            [vpara_max, vperp1_max, vperp2_max])
        vdf_fit_un_pc_tmp = 1.e-3 * dens_pc_lst[i_time] / (
                wthp_pc_lst[i_time] ** 2 * wthz_pc_lst[i_time] * np.pi ** 1.5) * \
                            np.exp(-(vpara - ux_tmp) ** 2 / wthz_pc_lst[i_time] ** 2) * \
                            np.exp(-(vperp1 - uy_tmp) ** 2 / wthp_pc_lst[i_time] ** 2) * \
                            np.exp(-(vperp2 - uz_tmp) ** 2 / wthp_pc_lst[i_time] ** 2) + \
                            1.e-3 * dens_pb_lst[i_time] / (
                                    wthp_pb_lst[i_time] ** 2 * wthz_pb_lst[i_time] * np.pi ** 1.5) * \
                            np.exp(-(vpara - ux_pb_lst[i_time]) ** 2 / wthz_pb_lst[i_time] ** 2) * \
                            np.exp(-(vperp1 - uy_pb_lst[i_time]) ** 2 / wthp_pb_lst[i_time] ** 2) * \
                            np.exp(-(vperp2 - uz_pb_lst[i_time]) ** 2 / wthp_pb_lst[i_time] ** 2)
        chi_un_pc_lst.append(np.nansum((vdf - vdf_fit_un_pc_tmp) ** 2 / sigma))
    for i_chi in range(num_var):
        ux_tmp, uy_tmp, uz_tmp = np.dot(rtn2mfa, np.vstack(
            [var_ur_pb_lst[i_chi], ut_pb_lst[i_time], un_pb_lst[i_time]])) - np.vstack(
            [vpara_max, vperp1_max, vperp2_max])
        vdf_fit_ur_pb_tmp = 1.e-3 * dens_pc_lst[i_time] / (
                wthp_pc_lst[i_time] ** 2 * wthz_pc_lst[i_time] * np.pi ** 1.5) * \
                            np.exp(-(vpara - ux_pc_lst[i_time]) ** 2 / wthz_pc_lst[i_time] ** 2) * \
                            np.exp(-(vperp1 - uy_pc_lst[i_time]) ** 2 / wthp_pc_lst[i_time] ** 2) * \
                            np.exp(-(vperp2 - uz_pc_lst[i_time]) ** 2 / wthp_pc_lst[i_time] ** 2) + \
                            1.e-3 * dens_pb_lst[i_time] / (
                                    wthp_pb_lst[i_time] ** 2 * wthz_pb_lst[i_time] * np.pi ** 1.5) * \
                            np.exp(-(vpara - ux_tmp) ** 2 / wthz_pb_lst[i_time] ** 2) * \
                            np.exp(-(vperp1 - uy_tmp) ** 2 / wthp_pb_lst[i_time] ** 2) * \
                            np.exp(-(vperp2 - uz_tmp) ** 2 / wthp_pb_lst[i_time] ** 2)
        chi_ur_pb_lst.append(np.nansum((vdf - vdf_fit_ur_pb_tmp) ** 2 / sigma))
    for i_chi in range(num_var):
        ux_tmp, uy_tmp, uz_tmp = np.dot(rtn2mfa, np.vstack(
            [ur_pb_lst[i_time], var_ut_pb_lst[i_chi], un_pb_lst[i_time]])) - np.vstack(
            [vpara_max, vperp1_max, vperp2_max])
        vdf_fit_ut_pb_tmp = 1.e-3 * dens_pc_lst[i_time] / (
                wthp_pc_lst[i_time] ** 2 * wthz_pc_lst[i_time] * np.pi ** 1.5) * \
                            np.exp(-(vpara - ux_pc_lst[i_time]) ** 2 / wthz_pc_lst[i_time] ** 2) * \
                            np.exp(-(vperp1 - uy_pc_lst[i_time]) ** 2 / wthp_pc_lst[i_time] ** 2) * \
                            np.exp(-(vperp2 - uz_pc_lst[i_time]) ** 2 / wthp_pc_lst[i_time] ** 2) + \
                            1.e-3 * dens_pb_lst[i_time] / (
                                    wthp_pb_lst[i_time] ** 2 * wthz_pb_lst[i_time] * np.pi ** 1.5) * \
                            np.exp(-(vpara - ux_tmp) ** 2 / wthz_pb_lst[i_time] ** 2) * \
                            np.exp(-(vperp1 - uy_tmp) ** 2 / wthp_pb_lst[i_time] ** 2) * \
                            np.exp(-(vperp2 - uz_tmp) ** 2 / wthp_pb_lst[i_time] ** 2)
        chi_ut_pb_lst.append(np.nansum((vdf - vdf_fit_ut_pb_tmp) ** 2 / sigma))
    for i_chi in range(num_var):
        ux_tmp, uy_tmp, uz_tmp = np.dot(rtn2mfa, np.vstack(
            [ur_pb_lst[i_time], ut_pb_lst[i_time], var_un_pb_lst[i_chi]])) - np.vstack(
            [vpara_max, vperp1_max, vperp2_max])
        vdf_fit_un_pb_tmp = 1.e-3 * dens_pc_lst[i_time] / (
                wthp_pc_lst[i_time] ** 2 * wthz_pc_lst[i_time] * np.pi ** 1.5) * \
                            np.exp(-(vpara - ux_pc_lst[i_time]) ** 2 / wthz_pc_lst[i_time] ** 2) * \
                            np.exp(-(vperp1 - uy_pc_lst[i_time]) ** 2 / wthp_pc_lst[i_time] ** 2) * \
                            np.exp(-(vperp2 - uz_pc_lst[i_time]) ** 2 / wthp_pc_lst[i_time] ** 2) + \
                            1.e-3 * dens_pb_lst[i_time] / (
                                    wthp_pb_lst[i_time] ** 2 * wthz_pb_lst[i_time] * np.pi ** 1.5) * \
                            np.exp(-(vpara - ux_tmp) ** 2 / wthz_pb_lst[i_time] ** 2) * \
                            np.exp(-(vperp1 - uy_tmp) ** 2 / wthp_pb_lst[i_time] ** 2) * \
                            np.exp(-(vperp2 - uz_tmp) ** 2 / wthp_pb_lst[i_time] ** 2)
        chi_un_pb_lst.append(np.nansum((vdf - vdf_fit_un_pb_tmp) ** 2 / sigma))
    if is_find_chi_square_min == 1:
        min_ur_pc_ind, min_ur_pb_ind = np.argmin(chi_ur_pc_lst), np.argmin(chi_ur_pb_lst)
        min_ut_pc_ind, min_ut_pb_ind = np.argmin(chi_ut_pc_lst), np.argmin(chi_ut_pb_lst)
        min_un_pc_ind, min_un_pb_ind = np.argmin(chi_un_pc_lst), np.argmin(chi_un_pb_lst)
        ur_pc_chi_min, ur_pb_chi_min = var_ur_pc_lst[min_ur_pc_ind], var_ur_pb_lst[min_ur_pb_ind]
        ut_pc_chi_min, ut_pb_chi_min = var_ut_pc_lst[min_ut_pc_ind], var_ut_pb_lst[min_ut_pb_ind]
        un_pc_chi_min, un_pb_chi_min = var_un_pc_lst[min_un_pc_ind], var_un_pb_lst[min_un_pb_ind]
        # recalculate X2 for individual parameters around minimum χ2 (rtn coordinates)
        chi_ur_pc_lst, chi_ur_pb_lst = [], []
        chi_ut_pc_lst, chi_ut_pb_lst = [], []
        chi_un_pc_lst, chi_un_pb_lst = [], []
        var_ur_pc_lst, var_ur_pb_lst = np.linspace(ur_pc_chi_min - 15, ur_pc_chi_min + 15, num_var), np.linspace(
            ur_pb_chi_min - 15, ur_pb_chi_min + 15, num_var)
        var_ut_pc_lst, var_ut_pb_lst = np.linspace(ut_pc_chi_min - 15, ut_pc_chi_min + 15, num_var), np.linspace(
            ut_pb_chi_min - 15, ut_pb_chi_min + 15, num_var)
        var_un_pc_lst, var_un_pb_lst = np.linspace(un_pc_chi_min - 15, un_pc_chi_min + 15, num_var), np.linspace(
            un_pb_chi_min - 15, un_pb_chi_min + 15, num_var)
        for i_chi in range(num_var):
            ux_tmp, uy_tmp, uz_tmp = np.dot(rtn2mfa, np.vstack(
                [var_ur_pc_lst[i_chi], ut_pc_lst[i_time], un_pc_lst[i_time]])) - np.vstack(
                [vpara_max, vperp1_max, vperp2_max])
            vdf_fit_ur_pc_tmp = 1.e-3 * dens_pc_lst[i_time] / (
                    wthp_pc_lst[i_time] ** 2 * wthz_pc_lst[i_time] * np.pi ** 1.5) * \
                                np.exp(-(vpara - ux_tmp) ** 2 / wthz_pc_lst[i_time] ** 2) * \
                                np.exp(-(vperp1 - uy_tmp) ** 2 / wthp_pc_lst[i_time] ** 2) * \
                                np.exp(-(vperp2 - uz_tmp) ** 2 / wthp_pc_lst[i_time] ** 2) + \
                                1.e-3 * dens_pb_lst[i_time] / (
                                        wthp_pb_lst[i_time] ** 2 * wthz_pb_lst[i_time] * np.pi ** 1.5) * \
                                np.exp(-(vpara - ux_pb_lst[i_time]) ** 2 / wthz_pb_lst[i_time] ** 2) * \
                                np.exp(-(vperp1 - uy_pb_lst[i_time]) ** 2 / wthp_pb_lst[i_time] ** 2) * \
                                np.exp(-(vperp2 - uz_pb_lst[i_time]) ** 2 / wthp_pb_lst[i_time] ** 2)
            chi_ur_pc_lst.append(np.nansum((vdf - vdf_fit_ur_pc_tmp) ** 2 / sigma))
        for i_chi in range(num_var):
            ux_tmp, uy_tmp, uz_tmp = np.dot(rtn2mfa, np.vstack(
                [ur_pc_lst[i_time], var_ut_pc_lst[i_chi], un_pc_lst[i_time]])) - np.vstack(
                [vpara_max, vperp1_max, vperp2_max])
            vdf_fit_ut_pc_tmp = 1.e-3 * dens_pc_lst[i_time] / (
                    wthp_pc_lst[i_time] ** 2 * wthz_pc_lst[i_time] * np.pi ** 1.5) * \
                                np.exp(-(vpara - ux_tmp) ** 2 / wthz_pc_lst[i_time] ** 2) * \
                                np.exp(-(vperp1 - uy_tmp) ** 2 / wthp_pc_lst[i_time] ** 2) * \
                                np.exp(-(vperp2 - uz_tmp) ** 2 / wthp_pc_lst[i_time] ** 2) + \
                                1.e-3 * dens_pb_lst[i_time] / (
                                        wthp_pb_lst[i_time] ** 2 * wthz_pb_lst[i_time] * np.pi ** 1.5) * \
                                np.exp(-(vpara - ux_pb_lst[i_time]) ** 2 / wthz_pb_lst[i_time] ** 2) * \
                                np.exp(-(vperp1 - uy_pb_lst[i_time]) ** 2 / wthp_pb_lst[i_time] ** 2) * \
                                np.exp(-(vperp2 - uz_pb_lst[i_time]) ** 2 / wthp_pb_lst[i_time] ** 2)
            chi_ut_pc_lst.append(np.nansum((vdf - vdf_fit_ut_pc_tmp) ** 2 / sigma))
        for i_chi in range(num_var):
            ux_tmp, uy_tmp, uz_tmp = np.dot(rtn2mfa, np.vstack(
                [ur_pc_lst[i_time], ut_pc_lst[i_time], var_un_pc_lst[i_chi]])) - np.vstack(
                [vpara_max, vperp1_max, vperp2_max])
            vdf_fit_un_pc_tmp = 1.e-3 * dens_pc_lst[i_time] / (
                    wthp_pc_lst[i_time] ** 2 * wthz_pc_lst[i_time] * np.pi ** 1.5) * \
                                np.exp(-(vpara - ux_tmp) ** 2 / wthz_pc_lst[i_time] ** 2) * \
                                np.exp(-(vperp1 - uy_tmp) ** 2 / wthp_pc_lst[i_time] ** 2) * \
                                np.exp(-(vperp2 - uz_tmp) ** 2 / wthp_pc_lst[i_time] ** 2) + \
                                1.e-3 * dens_pb_lst[i_time] / (
                                        wthp_pb_lst[i_time] ** 2 * wthz_pb_lst[i_time] * np.pi ** 1.5) * \
                                np.exp(-(vpara - ux_pb_lst[i_time]) ** 2 / wthz_pb_lst[i_time] ** 2) * \
                                np.exp(-(vperp1 - uy_pb_lst[i_time]) ** 2 / wthp_pb_lst[i_time] ** 2) * \
                                np.exp(-(vperp2 - uz_pb_lst[i_time]) ** 2 / wthp_pb_lst[i_time] ** 2)
            chi_un_pc_lst.append(np.nansum((vdf - vdf_fit_un_pc_tmp) ** 2 / sigma))
        for i_chi in range(num_var):
            ux_tmp, uy_tmp, uz_tmp = np.dot(rtn2mfa, np.vstack(
                [var_ur_pb_lst[i_chi], ut_pb_lst[i_time], un_pb_lst[i_time]])) - np.vstack(
                [vpara_max, vperp1_max, vperp2_max])
            vdf_fit_ur_pb_tmp = 1.e-3 * dens_pc_lst[i_time] / (
                    wthp_pc_lst[i_time] ** 2 * wthz_pc_lst[i_time] * np.pi ** 1.5) * \
                                np.exp(-(vpara - ux_pc_lst[i_time]) ** 2 / wthz_pc_lst[i_time] ** 2) * \
                                np.exp(-(vperp1 - uy_pc_lst[i_time]) ** 2 / wthp_pc_lst[i_time] ** 2) * \
                                np.exp(-(vperp2 - uz_pc_lst[i_time]) ** 2 / wthp_pc_lst[i_time] ** 2) + \
                                1.e-3 * dens_pb_lst[i_time] / (
                                        wthp_pb_lst[i_time] ** 2 * wthz_pb_lst[i_time] * np.pi ** 1.5) * \
                                np.exp(-(vpara - ux_tmp) ** 2 / wthz_pb_lst[i_time] ** 2) * \
                                np.exp(-(vperp1 - uy_tmp) ** 2 / wthp_pb_lst[i_time] ** 2) * \
                                np.exp(-(vperp2 - uz_tmp) ** 2 / wthp_pb_lst[i_time] ** 2)
            chi_ur_pb_lst.append(np.nansum((vdf - vdf_fit_ur_pb_tmp) ** 2 / sigma))
        for i_chi in range(num_var):
            ux_tmp, uy_tmp, uz_tmp = np.dot(rtn2mfa, np.vstack(
                [ur_pb_lst[i_time], var_ut_pb_lst[i_chi], un_pb_lst[i_time]])) - np.vstack(
                [vpara_max, vperp1_max, vperp2_max])
            vdf_fit_ut_pb_tmp = 1.e-3 * dens_pc_lst[i_time] / (
                    wthp_pc_lst[i_time] ** 2 * wthz_pc_lst[i_time] * np.pi ** 1.5) * \
                                np.exp(-(vpara - ux_pc_lst[i_time]) ** 2 / wthz_pc_lst[i_time] ** 2) * \
                                np.exp(-(vperp1 - uy_pc_lst[i_time]) ** 2 / wthp_pc_lst[i_time] ** 2) * \
                                np.exp(-(vperp2 - uz_pc_lst[i_time]) ** 2 / wthp_pc_lst[i_time] ** 2) + \
                                1.e-3 * dens_pb_lst[i_time] / (
                                        wthp_pb_lst[i_time] ** 2 * wthz_pb_lst[i_time] * np.pi ** 1.5) * \
                                np.exp(-(vpara - ux_tmp) ** 2 / wthz_pb_lst[i_time] ** 2) * \
                                np.exp(-(vperp1 - uy_tmp) ** 2 / wthp_pb_lst[i_time] ** 2) * \
                                np.exp(-(vperp2 - uz_tmp) ** 2 / wthp_pb_lst[i_time] ** 2)
            chi_ut_pb_lst.append(np.nansum((vdf - vdf_fit_ut_pb_tmp) ** 2 / sigma))
        for i_chi in range(num_var):
            ux_tmp, uy_tmp, uz_tmp = np.dot(rtn2mfa, np.vstack(
                [ur_pb_lst[i_time], ut_pb_lst[i_time], var_un_pb_lst[i_chi]])) - np.vstack(
                [vpara_max, vperp1_max, vperp2_max])
            vdf_fit_un_pb_tmp = 1.e-3 * dens_pc_lst[i_time] / (
                    wthp_pc_lst[i_time] ** 2 * wthz_pc_lst[i_time] * np.pi ** 1.5) * \
                                np.exp(-(vpara - ux_pc_lst[i_time]) ** 2 / wthz_pc_lst[i_time] ** 2) * \
                                np.exp(-(vperp1 - uy_pc_lst[i_time]) ** 2 / wthp_pc_lst[i_time] ** 2) * \
                                np.exp(-(vperp2 - uz_pc_lst[i_time]) ** 2 / wthp_pc_lst[i_time] ** 2) + \
                                1.e-3 * dens_pb_lst[i_time] / (
                                        wthp_pb_lst[i_time] ** 2 * wthz_pb_lst[i_time] * np.pi ** 1.5) * \
                                np.exp(-(vpara - ux_tmp) ** 2 / wthz_pb_lst[i_time] ** 2) * \
                                np.exp(-(vperp1 - uy_tmp) ** 2 / wthp_pb_lst[i_time] ** 2) * \
                                np.exp(-(vperp2 - uz_tmp) ** 2 / wthp_pb_lst[i_time] ** 2)
            chi_un_pb_lst.append(np.nansum((vdf - vdf_fit_un_pb_tmp) ** 2 / sigma))
    # determine the uncertainties of the parameters using (χ2 - Δχ2, χ2 + Δχ2) [confidence level: 99.99$]
    # according to Numerical Recipes The Art of Scientific Computing
    delta_chi2 = 6.63
    # bulk velocity (u_r) proton core
    chi_fit = np.nansum((vdf - vdf_fit - mean_times_ori_sd) ** 2 / sigma)
    ind = np.argmin(np.abs(chi_ur_pc_lst - chi_fit))
    min_chi, max_chi = chi_ur_pc_lst[ind] - delta_chi2, chi_ur_pc_lst[ind] + delta_chi2
    ind_with = np.where(np.logical_and(chi_ur_pc_lst > min_chi, chi_ur_pc_lst < max_chi))
    ur_pc_low_lst.append(var_ur_pc_lst[np.min(ind_with)])
    ur_pc_up_lst.append(var_ur_pc_lst[np.max(ind_with)])
    delta_ur_pc_low_lst.append(ur_pc_chi_min - var_ur_pc_lst[np.min(ind_with)])
    delta_ur_pc_up_lst.append(var_ur_pc_lst[np.max(ind_with)] - ur_pc_chi_min)
    # bulk velocity (u_t) proton core
    chi_fit = np.nansum((vdf - vdf_fit - mean_times_ori_sd) ** 2 / sigma)
    ind = np.argmin(np.abs(chi_ut_pc_lst - chi_fit))
    min_chi, max_chi = chi_ut_pc_lst[ind] - delta_chi2, chi_ut_pc_lst[ind] + delta_chi2
    ind_with = np.where(np.logical_and(chi_ut_pc_lst > min_chi, chi_ut_pc_lst < max_chi))
    ut_pc_low_lst.append(var_ut_pc_lst[np.min(ind_with)])
    ut_pc_up_lst.append(var_ut_pc_lst[np.max(ind_with)])
    delta_ut_pc_low_lst.append(ut_pc_chi_min - var_ut_pc_lst[np.min(ind_with)])
    delta_ut_pc_up_lst.append(var_ut_pc_lst[np.max(ind_with)] - ut_pc_chi_min)
    # bulk velocity (u_n) proton core
    chi_fit = np.nansum((vdf - vdf_fit - mean_times_ori_sd) ** 2 / sigma)
    ind = np.argmin(np.abs(chi_un_pc_lst - chi_fit))
    min_chi, max_chi = chi_un_pc_lst[ind] - delta_chi2, chi_un_pc_lst[ind] + delta_chi2
    ind_with = np.where(np.logical_and(chi_un_pc_lst > min_chi, chi_un_pc_lst < max_chi))
    un_pc_low_lst.append(var_un_pc_lst[np.min(ind_with)])
    un_pc_up_lst.append(var_un_pc_lst[np.max(ind_with)])
    delta_un_pc_low_lst.append(un_pc_chi_min - var_un_pc_lst[np.min(ind_with)])
    delta_un_pc_up_lst.append(var_un_pc_lst[np.max(ind_with)] - un_pc_chi_min)
    # bulk velocity (u_r) proton beam
    chi_fit = np.nansum((vdf - vdf_fit - mean_times_ori_sd) ** 2 / sigma)
    ind = np.argmin(np.abs(chi_ur_pb_lst - chi_fit))
    min_chi, max_chi = chi_ur_pb_lst[ind] - delta_chi2, chi_ur_pb_lst[ind] + delta_chi2
    ind_with = np.where(np.logical_and(chi_ur_pb_lst > min_chi, chi_ur_pb_lst < max_chi))
    ur_pb_low_lst.append(var_ur_pb_lst[np.min(ind_with)])
    ur_pb_up_lst.append(var_ur_pb_lst[np.max(ind_with)])
    delta_ur_pb_low_lst.append(ur_pb_chi_min - var_ur_pb_lst[np.min(ind_with)])
    delta_ur_pb_up_lst.append(var_ur_pb_lst[np.max(ind_with)] - ur_pb_chi_min)
    # bulk velocity (u_t) proton beam
    chi_fit = np.nansum((vdf - vdf_fit - mean_times_ori_sd) ** 2 / sigma)
    ind = np.argmin(np.abs(chi_ut_pb_lst - chi_fit))
    min_chi, max_chi = chi_ut_pb_lst[ind] - delta_chi2, chi_ut_pb_lst[ind] + delta_chi2
    ind_with = np.where(np.logical_and(chi_ut_pb_lst > min_chi, chi_ut_pb_lst < max_chi))
    ut_pb_low_lst.append(var_ut_pb_lst[np.min(ind_with)])
    ut_pb_up_lst.append(var_ut_pb_lst[np.max(ind_with)])
    delta_ut_pb_low_lst.append(ut_pb_chi_min - var_ut_pb_lst[np.min(ind_with)])
    delta_ut_pb_up_lst.append(var_ut_pb_lst[np.max(ind_with)] - ut_pb_chi_min)
    # bulk velocity (u_n) proton beam
    chi_fit = np.nansum((vdf - vdf_fit - mean_times_ori_sd) ** 2 / sigma)
    ind = np.argmin(np.abs(chi_un_pb_lst - chi_fit))
    min_chi, max_chi = chi_un_pb_lst[ind] - delta_chi2, chi_un_pb_lst[ind] + delta_chi2
    ind_with = np.where(np.logical_and(chi_un_pb_lst > min_chi, chi_un_pb_lst < max_chi))
    un_pb_low_lst.append(var_un_pb_lst[np.min(ind_with)])
    un_pb_up_lst.append(var_un_pb_lst[np.max(ind_with)])
    delta_un_pb_low_lst.append(un_pb_chi_min - var_un_pb_lst[np.min(ind_with)])
    delta_un_pb_up_lst.append(var_un_pb_lst[np.max(ind_with)] - un_pb_chi_min)
    # calculate Δχ2 of individual parameter with other parameters fixed in mfa coordinates
    for i_chi in range(num_var):
        vdf_fit_n_pc_tmp = 1.e-3 * var_n_pc_lst[i_chi] / (
                wthp_pc_lst[i_time] ** 2 * wthz_pc_lst[i_time] * np.pi ** 1.5) * \
                           np.exp(-(vpara - ux_pc_lst[i_time]) ** 2 / wthz_pc_lst[i_time] ** 2) * \
                           np.exp(-(vperp1 - uy_pc_lst[i_time]) ** 2 / wthp_pc_lst[i_time] ** 2) * \
                           np.exp(-(vperp2 - uz_pc_lst[i_time]) ** 2 / wthp_pc_lst[i_time] ** 2) + \
                           1.e-3 * dens_pb_lst[i_time] / (
                                   wthp_pb_lst[i_time] ** 2 * wthz_pb_lst[i_time] * np.pi ** 1.5) * \
                           np.exp(-(vpara - ux_pb_lst[i_time]) ** 2 / wthz_pb_lst[i_time] ** 2) * \
                           np.exp(-(vperp1 - uy_pb_lst[i_time]) ** 2 / wthp_pb_lst[i_time] ** 2) * \
                           np.exp(-(vperp2 - uz_pb_lst[i_time]) ** 2 / wthp_pb_lst[i_time] ** 2)
        vdf_fit_n_pb_tmp = 1.e-3 * dens_pc_lst[i_time] / (
                wthp_pc_lst[i_time] ** 2 * wthz_pc_lst[i_time] * np.pi ** 1.5) * \
                           np.exp(-(vpara - ux_pc_lst[i_time]) ** 2 / wthz_pc_lst[i_time] ** 2) * \
                           np.exp(-(vperp1 - uy_pc_lst[i_time]) ** 2 / wthp_pc_lst[i_time] ** 2) * \
                           np.exp(-(vperp2 - uz_pc_lst[i_time]) ** 2 / wthp_pc_lst[i_time] ** 2) + \
                           1.e-3 * var_n_pb_lst[i_chi] / (
                                   wthp_pb_lst[i_time] ** 2 * wthz_pb_lst[i_time] * np.pi ** 1.5) * \
                           np.exp(-(vpara - ux_pb_lst[i_time]) ** 2 / wthz_pb_lst[i_time] ** 2) * \
                           np.exp(-(vperp1 - uy_pb_lst[i_time]) ** 2 / wthp_pb_lst[i_time] ** 2) * \
                           np.exp(-(vperp2 - uz_pb_lst[i_time]) ** 2 / wthp_pb_lst[i_time] ** 2)
        vdf_fit_wthp_pc_tmp = 1.e-3 * dens_pc_lst[i_time] / (
                var_wthp_pc_lst[i_chi] ** 2 * wthz_pc_lst[i_time] * np.pi ** 1.5) * \
                              np.exp(-(vpara - ux_pc_lst[i_time]) ** 2 / wthz_pc_lst[i_time] ** 2) * \
                              np.exp(-(vperp1 - uy_pc_lst[i_time]) ** 2 / var_wthp_pc_lst[i_chi] ** 2) * \
                              np.exp(-(vperp2 - uz_pc_lst[i_time]) ** 2 / var_wthp_pc_lst[i_chi] ** 2) + \
                              1.e-3 * dens_pb_lst[i_time] / (
                                      wthp_pb_lst[i_time] ** 2 * wthz_pb_lst[i_time] * np.pi ** 1.5) * \
                              np.exp(-(vpara - ux_pb_lst[i_time]) ** 2 / wthz_pb_lst[i_time] ** 2) * \
                              np.exp(-(vperp1 - uy_pb_lst[i_time]) ** 2 / wthp_pb_lst[i_time] ** 2) * \
                              np.exp(-(vperp2 - uz_pb_lst[i_time]) ** 2 / wthp_pb_lst[i_time] ** 2)
        vdf_fit_wthp_pb_tmp = 1.e-3 * dens_pc_lst[i_time] / (
                wthp_pc_lst[i_time] ** 2 * wthz_pc_lst[i_time] * np.pi ** 1.5) * \
                              np.exp(-(vpara - ux_pc_lst[i_time]) ** 2 / wthz_pc_lst[i_time] ** 2) * \
                              np.exp(-(vperp1 - uy_pc_lst[i_time]) ** 2 / wthp_pc_lst[i_time] ** 2) * \
                              np.exp(-(vperp2 - uz_pc_lst[i_time]) ** 2 / wthp_pc_lst[i_time] ** 2) + \
                              1.e-3 * dens_pb_lst[i_time] / (
                                      var_wthp_pc_lst[i_chi] ** 2 * wthz_pb_lst[i_time] * np.pi ** 1.5) * \
                              np.exp(-(vpara - ux_pb_lst[i_time]) ** 2 / wthz_pb_lst[i_time] ** 2) * \
                              np.exp(-(vperp1 - uy_pb_lst[i_time]) ** 2 / var_wthp_pc_lst[i_chi] ** 2) * \
                              np.exp(-(vperp2 - uz_pb_lst[i_time]) ** 2 / var_wthp_pc_lst[i_chi] ** 2)
        vdf_fit_wthz_pc_tmp = 1.e-3 * dens_pc_lst[i_time] / (
                wthp_pc_lst[i_time] ** 2 * var_wthz_pc_lst[i_chi] * np.pi ** 1.5) * \
                              np.exp(-(vpara - ux_pc_lst[i_time]) ** 2 / var_wthz_pc_lst[i_chi] ** 2) * \
                              np.exp(-(vperp1 - uy_pc_lst[i_time]) ** 2 / wthp_pc_lst[i_time] ** 2) * \
                              np.exp(-(vperp2 - uz_pc_lst[i_time]) ** 2 / wthp_pc_lst[i_time] ** 2) + \
                              1.e-3 * dens_pb_lst[i_time] / (
                                      wthp_pb_lst[i_time] ** 2 * wthz_pb_lst[i_time] * np.pi ** 1.5) * \
                              np.exp(-(vpara - ux_pb_lst[i_time]) ** 2 / wthz_pb_lst[i_time] ** 2) * \
                              np.exp(-(vperp1 - uy_pb_lst[i_time]) ** 2 / wthp_pb_lst[i_time] ** 2) * \
                              np.exp(-(vperp2 - uz_pb_lst[i_time]) ** 2 / wthp_pb_lst[i_time] ** 2)
        vdf_fit_wthz_pb_tmp = 1.e-3 * dens_pc_lst[i_time] / (
                wthp_pc_lst[i_time] ** 2 * wthz_pc_lst[i_time] * np.pi ** 1.5) * \
                              np.exp(-(vpara - ux_pc_lst[i_time]) ** 2 / wthz_pc_lst[i_time] ** 2) * \
                              np.exp(-(vperp1 - uy_pc_lst[i_time]) ** 2 / wthp_pc_lst[i_time] ** 2) * \
                              np.exp(-(vperp2 - uz_pc_lst[i_time]) ** 2 / wthp_pc_lst[i_time] ** 2) + \
                              1.e-3 * dens_pb_lst[i_time] / (
                                      wthp_pb_lst[i_time] ** 2 * var_wthz_pb_lst[i_chi] * np.pi ** 1.5) * \
                              np.exp(-(vpara - ux_pb_lst[i_time]) ** 2 / var_wthz_pb_lst[i_chi] ** 2) * \
                              np.exp(-(vperp1 - uy_pb_lst[i_time]) ** 2 / wthp_pb_lst[i_time] ** 2) * \
                              np.exp(-(vperp2 - uz_pb_lst[i_time]) ** 2 / wthp_pb_lst[i_time] ** 2)
        vdf_fit_ux_pc_tmp = 1.e-3 * dens_pc_lst[i_time] / (
                wthp_pc_lst[i_time] ** 2 * wthz_pc_lst[i_time] * np.pi ** 1.5) * \
                            np.exp(-(vpara - var_ux_pc_lst[i_chi]) ** 2 / wthz_pc_lst[i_time] ** 2) * \
                            np.exp(-(vperp1 - uy_pc_lst[i_time]) ** 2 / wthp_pc_lst[i_time] ** 2) * \
                            np.exp(-(vperp2 - uz_pc_lst[i_time]) ** 2 / wthp_pc_lst[i_time] ** 2) + \
                            1.e-3 * dens_pb_lst[i_time] / (
                                    wthp_pb_lst[i_time] ** 2 * wthz_pb_lst[i_time] * np.pi ** 1.5) * \
                            np.exp(-(vpara - ux_pb_lst[i_time]) ** 2 / wthz_pb_lst[i_time] ** 2) * \
                            np.exp(-(vperp1 - uy_pb_lst[i_time]) ** 2 / wthp_pb_lst[i_time] ** 2) * \
                            np.exp(-(vperp2 - uz_pb_lst[i_time]) ** 2 / wthp_pb_lst[i_time] ** 2)
        vdf_fit_ux_pb_tmp = 1.e-3 * dens_pc_lst[i_time] / (
                wthp_pc_lst[i_time] ** 2 * wthz_pc_lst[i_time] * np.pi ** 1.5) * \
                            np.exp(-(vpara - ux_pc_lst[i_time]) ** 2 / wthz_pc_lst[i_time] ** 2) * \
                            np.exp(-(vperp1 - uy_pc_lst[i_time]) ** 2 / wthp_pc_lst[i_time] ** 2) * \
                            np.exp(-(vperp2 - uz_pc_lst[i_time]) ** 2 / wthp_pc_lst[i_time] ** 2) + \
                            1.e-3 * dens_pb_lst[i_time] / (
                                    wthp_pb_lst[i_time] ** 2 * wthz_pb_lst[i_time] * np.pi ** 1.5) * \
                            np.exp(-(vpara - var_ux_pb_lst[i_chi]) ** 2 / wthz_pb_lst[i_time] ** 2) * \
                            np.exp(-(vperp1 - uy_pb_lst[i_time]) ** 2 / wthp_pb_lst[i_time] ** 2) * \
                            np.exp(-(vperp2 - uz_pb_lst[i_time]) ** 2 / wthp_pb_lst[i_time] ** 2)
        vdf_fit_uy_pc_tmp = 1.e-3 * dens_pc_lst[i_time] / (
                wthp_pc_lst[i_time] ** 2 * wthz_pc_lst[i_time] * np.pi ** 1.5) * \
                            np.exp(-(vpara - ux_pc_lst[i_time]) ** 2 / wthz_pc_lst[i_time] ** 2) * \
                            np.exp(-(vperp1 - var_uy_pc_lst[i_chi]) ** 2 / wthp_pc_lst[i_time] ** 2) * \
                            np.exp(-(vperp2 - uz_pc_lst[i_time]) ** 2 / wthp_pc_lst[i_time] ** 2) + \
                            1.e-3 * dens_pb_lst[i_time] / (
                                    wthp_pb_lst[i_time] ** 2 * wthz_pb_lst[i_time] * np.pi ** 1.5) * \
                            np.exp(-(vpara - ux_pb_lst[i_time]) ** 2 / wthz_pb_lst[i_time] ** 2) * \
                            np.exp(-(vperp1 - uy_pb_lst[i_time]) ** 2 / wthp_pb_lst[i_time] ** 2) * \
                            np.exp(-(vperp2 - uz_pb_lst[i_time]) ** 2 / wthp_pb_lst[i_time] ** 2)
        vdf_fit_uy_pb_tmp = 1.e-3 * dens_pc_lst[i_time] / (
                wthp_pc_lst[i_time] ** 2 * wthz_pc_lst[i_time] * np.pi ** 1.5) * \
                            np.exp(-(vpara - ux_pc_lst[i_time]) ** 2 / wthz_pc_lst[i_time] ** 2) * \
                            np.exp(-(vperp1 - uy_pc_lst[i_time]) ** 2 / wthp_pc_lst[i_time] ** 2) * \
                            np.exp(-(vperp2 - uz_pc_lst[i_time]) ** 2 / wthp_pc_lst[i_time] ** 2) + \
                            1.e-3 * dens_pb_lst[i_time] / (
                                    wthp_pb_lst[i_time] ** 2 * wthz_pb_lst[i_time] * np.pi ** 1.5) * \
                            np.exp(-(vpara - ux_pb_lst[i_time]) ** 2 / wthz_pb_lst[i_time] ** 2) * \
                            np.exp(-(vperp1 - var_uy_pb_lst[i_chi]) ** 2 / wthp_pb_lst[i_time] ** 2) * \
                            np.exp(-(vperp2 - uz_pb_lst[i_time]) ** 2 / wthp_pb_lst[i_time] ** 2)
        vdf_fit_uz_pc_tmp = 1.e-3 * dens_pc_lst[i_time] / (
                wthp_pc_lst[i_time] ** 2 * wthz_pc_lst[i_time] * np.pi ** 1.5) * \
                            np.exp(-(vpara - ux_pc_lst[i_time]) ** 2 / wthz_pc_lst[i_time] ** 2) * \
                            np.exp(-(vperp1 - uy_pc_lst[i_time]) ** 2 / wthp_pc_lst[i_time] ** 2) * \
                            np.exp(-(vperp2 - var_uz_pc_lst[i_chi]) ** 2 / wthp_pc_lst[i_time] ** 2) + \
                            1.e-3 * dens_pb_lst[i_time] / (
                                    wthp_pb_lst[i_time] ** 2 * wthz_pb_lst[i_time] * np.pi ** 1.5) * \
                            np.exp(-(vpara - ux_pb_lst[i_time]) ** 2 / wthz_pb_lst[i_time] ** 2) * \
                            np.exp(-(vperp1 - uy_pb_lst[i_time]) ** 2 / wthp_pb_lst[i_time] ** 2) * \
                            np.exp(-(vperp2 - uz_pb_lst[i_time]) ** 2 / wthp_pb_lst[i_time] ** 2)
        vdf_fit_uz_pb_tmp = 1.e-3 * dens_pc_lst[i_time] / (
                wthp_pc_lst[i_time] ** 2 * wthz_pc_lst[i_time] * np.pi ** 1.5) * \
                            np.exp(-(vpara - ux_pc_lst[i_time]) ** 2 / wthz_pc_lst[i_time] ** 2) * \
                            np.exp(-(vperp1 - uy_pc_lst[i_time]) ** 2 / wthp_pc_lst[i_time] ** 2) * \
                            np.exp(-(vperp2 - uz_pc_lst[i_time]) ** 2 / wthp_pc_lst[i_time] ** 2) + \
                            1.e-3 * dens_pb_lst[i_time] / (
                                    wthp_pb_lst[i_time] ** 2 * wthz_pb_lst[i_time] * np.pi ** 1.5) * \
                            np.exp(-(vpara - ux_pb_lst[i_time]) ** 2 / wthz_pb_lst[i_time] ** 2) * \
                            np.exp(-(vperp1 - uy_pb_lst[i_time]) ** 2 / wthp_pb_lst[i_time] ** 2) * \
                            np.exp(-(vperp2 - var_uz_pb_lst[i_chi]) ** 2 / wthp_pb_lst[i_time] ** 2)
        chi_n_pc_lst.append(np.nansum((vdf - vdf_fit_n_pc_tmp - mean_times_ori_sd) ** 2 / sigma))
        chi_n_pb_lst.append(np.nansum((vdf - vdf_fit_n_pb_tmp - mean_times_ori_sd) ** 2 / sigma))
        chi_wthp_pc_lst.append(np.nansum((vdf - vdf_fit_wthp_pc_tmp - mean_times_ori_sd) ** 2 / sigma))
        chi_wthp_pb_lst.append(np.nansum((vdf - vdf_fit_wthp_pb_tmp - mean_times_ori_sd) ** 2 / sigma))
        chi_wthz_pc_lst.append(np.nansum((vdf - vdf_fit_wthz_pc_tmp - mean_times_ori_sd) ** 2 / sigma))
        chi_wthz_pb_lst.append(np.nansum((vdf - vdf_fit_wthz_pb_tmp - mean_times_ori_sd) ** 2 / sigma))
        chi_ux_pc_lst.append(np.nansum((vdf - vdf_fit_ux_pc_tmp - mean_times_ori_sd) ** 2 / sigma))
        chi_ux_pb_lst.append(np.nansum((vdf - vdf_fit_ux_pb_tmp - mean_times_ori_sd) ** 2 / sigma))
        chi_uy_pc_lst.append(np.nansum((vdf - vdf_fit_uy_pc_tmp - mean_times_ori_sd) ** 2 / sigma))
        chi_uy_pb_lst.append(np.nansum((vdf - vdf_fit_uy_pb_tmp - mean_times_ori_sd) ** 2 / sigma))
        chi_uz_pc_lst.append(np.nansum((vdf - vdf_fit_uz_pc_tmp - mean_times_ori_sd) ** 2 / sigma))
        chi_uz_pb_lst.append(np.nansum((vdf - vdf_fit_uz_pb_tmp - mean_times_ori_sd) ** 2 / sigma))
    if is_find_chi_square_min == 1:
        min_n_pc_ind, min_n_pb_ind = np.argmin(chi_n_pc_lst), np.argmin(chi_n_pb_lst)
        min_wthp_pc_ind, min_wthp_pb_ind = np.argmin(chi_wthp_pc_lst), np.argmin(chi_wthp_pb_lst)
        min_wthz_pc_ind, min_wthz_pb_ind = np.argmin(chi_wthz_pc_lst), np.argmin(chi_wthz_pb_lst)
        min_ux_pc_ind, min_ux_pb_ind = np.argmin(chi_ux_pc_lst), np.argmin(chi_ux_pb_lst)
        min_uy_pc_ind, min_uy_pb_ind = np.argmin(chi_uy_pc_lst), np.argmin(chi_uy_pb_lst)
        min_uz_pc_ind, min_uz_pb_ind = np.argmin(chi_uz_pc_lst), np.argmin(chi_uz_pb_lst)
        n_pc_chi_min, n_pb_chi_min = var_n_pc_lst[min_n_pc_ind], var_n_pb_lst[min_n_pb_ind]
        wthp_pc_chi_min, wthp_pb_chi_min = var_wthp_pc_lst[min_wthp_pc_ind], var_wthp_pb_lst[min_wthp_pb_ind]
        wthz_pc_chi_min, wthz_pb_chi_min = var_wthz_pc_lst[min_wthz_pc_ind], var_wthz_pb_lst[min_wthz_pb_ind]
        ux_pc_chi_min, ux_pb_chi_min = var_ux_pc_lst[min_ux_pc_ind], var_ux_pb_lst[min_ux_pb_ind]
        uy_pc_chi_min, uy_pb_chi_min = var_uy_pc_lst[min_uy_pc_ind], var_uy_pb_lst[min_uy_pb_ind]
        uz_pc_chi_min, uz_pb_chi_min = var_uz_pc_lst[min_uz_pc_ind], var_uz_pb_lst[min_uz_pb_ind]
        # recalculate X2 for individual parameters around minimum χ2 (rtn coordinates)
        chi_n_pc_lst, chi_n_pb_lst = [], []
        chi_wthp_pc_lst, chi_wthp_pb_lst, chi_wthz_pc_lst, chi_wthz_pb_lst = [], [], [], []
        chi_ux_pc_lst, chi_ux_pb_lst, chi_uy_pc_lst, chi_uy_pb_lst, chi_uz_pc_lst, chi_uz_pb_lst = [], [], [], [], [], []
        var_n_pc_lst = np.linspace(n_pc_chi_min - 2, n_pc_chi_min + 2, num_var)
        var_n_pb_lst = np.linspace(n_pb_chi_min - 4, n_pb_chi_min + 4, num_var)
        var_wthp_pc_lst = np.linspace(wthp_pc_chi_min - 15, wthp_pc_chi_min + 15, num_var)
        var_wthp_pb_lst = np.linspace(wthp_pb_chi_min - 35, wthp_pb_chi_min + 35, num_var)
        var_wthz_pc_lst = np.linspace(wthz_pc_chi_min - 15, wthz_pc_chi_min + 15, num_var)
        var_wthz_pb_lst = np.linspace(wthz_pb_chi_min - 35, wthz_pb_chi_min + 35, num_var)
        var_ux_pc_lst = np.linspace(ux_pc_chi_min - 15, ux_pc_chi_min + 15, num_var)
        var_ux_pb_lst = np.linspace(ux_pb_chi_min - 15, ux_pb_chi_min + 15, num_var)
        var_uy_pc_lst = np.linspace(uy_pc_chi_min - 15, uy_pc_chi_min + 15, num_var)
        var_uy_pb_lst = np.linspace(uy_pb_chi_min - 15, uy_pb_chi_min + 15, num_var)
        var_uz_pc_lst = np.linspace(uz_pc_chi_min - 15, uz_pc_chi_min + 15, num_var)
        var_uz_pb_lst = np.linspace(uz_pb_chi_min - 15, uz_pb_chi_min + 15, num_var)
        for i_chi in range(num_var):
            vdf_fit_n_pc_tmp = 1.e-3 * var_n_pc_lst[i_chi] / (
                    wthp_pc_lst[i_time] ** 2 * wthz_pc_lst[i_time] * np.pi ** 1.5) * \
                               np.exp(-(vpara - ux_pc_lst[i_time]) ** 2 / wthz_pc_lst[i_time] ** 2) * \
                               np.exp(-(vperp1 - uy_pc_lst[i_time]) ** 2 / wthp_pc_lst[i_time] ** 2) * \
                               np.exp(-(vperp2 - uz_pc_lst[i_time]) ** 2 / wthp_pc_lst[i_time] ** 2) + \
                               1.e-3 * dens_pb_lst[i_time] / (
                                       wthp_pb_lst[i_time] ** 2 * wthz_pb_lst[i_time] * np.pi ** 1.5) * \
                               np.exp(-(vpara - ux_pb_lst[i_time]) ** 2 / wthz_pb_lst[i_time] ** 2) * \
                               np.exp(-(vperp1 - uy_pb_lst[i_time]) ** 2 / wthp_pb_lst[i_time] ** 2) * \
                               np.exp(-(vperp2 - uz_pb_lst[i_time]) ** 2 / wthp_pb_lst[i_time] ** 2)
            vdf_fit_n_pb_tmp = 1.e-3 * dens_pc_lst[i_time] / (
                    wthp_pc_lst[i_time] ** 2 * wthz_pc_lst[i_time] * np.pi ** 1.5) * \
                               np.exp(-(vpara - ux_pc_lst[i_time]) ** 2 / wthz_pc_lst[i_time] ** 2) * \
                               np.exp(-(vperp1 - uy_pc_lst[i_time]) ** 2 / wthp_pc_lst[i_time] ** 2) * \
                               np.exp(-(vperp2 - uz_pc_lst[i_time]) ** 2 / wthp_pc_lst[i_time] ** 2) + \
                               1.e-3 * var_n_pb_lst[i_chi] / (
                                       wthp_pb_lst[i_time] ** 2 * wthz_pb_lst[i_time] * np.pi ** 1.5) * \
                               np.exp(-(vpara - ux_pb_lst[i_time]) ** 2 / wthz_pb_lst[i_time] ** 2) * \
                               np.exp(-(vperp1 - uy_pb_lst[i_time]) ** 2 / wthp_pb_lst[i_time] ** 2) * \
                               np.exp(-(vperp2 - uz_pb_lst[i_time]) ** 2 / wthp_pb_lst[i_time] ** 2)
            vdf_fit_wthp_pc_tmp = 1.e-3 * dens_pc_lst[i_time] / (
                    var_wthp_pc_lst[i_chi] ** 2 * wthz_pc_lst[i_time] * np.pi ** 1.5) * \
                                  np.exp(-(vpara - ux_pc_lst[i_time]) ** 2 / wthz_pc_lst[i_time] ** 2) * \
                                  np.exp(-(vperp1 - uy_pc_lst[i_time]) ** 2 / var_wthp_pc_lst[i_chi] ** 2) * \
                                  np.exp(-(vperp2 - uz_pc_lst[i_time]) ** 2 / var_wthp_pc_lst[i_chi] ** 2) + \
                                  1.e-3 * dens_pb_lst[i_time] / (
                                          wthp_pb_lst[i_time] ** 2 * wthz_pb_lst[i_time] * np.pi ** 1.5) * \
                                  np.exp(-(vpara - ux_pb_lst[i_time]) ** 2 / wthz_pb_lst[i_time] ** 2) * \
                                  np.exp(-(vperp1 - uy_pb_lst[i_time]) ** 2 / wthp_pb_lst[i_time] ** 2) * \
                                  np.exp(-(vperp2 - uz_pb_lst[i_time]) ** 2 / wthp_pb_lst[i_time] ** 2)
            vdf_fit_wthp_pb_tmp = 1.e-3 * dens_pc_lst[i_time] / (
                    wthp_pc_lst[i_time] ** 2 * wthz_pc_lst[i_time] * np.pi ** 1.5) * \
                                  np.exp(-(vpara - ux_pc_lst[i_time]) ** 2 / wthz_pc_lst[i_time] ** 2) * \
                                  np.exp(-(vperp1 - uy_pc_lst[i_time]) ** 2 / wthp_pc_lst[i_time] ** 2) * \
                                  np.exp(-(vperp2 - uz_pc_lst[i_time]) ** 2 / wthp_pc_lst[i_time] ** 2) + \
                                  1.e-3 * dens_pb_lst[i_time] / (
                                          var_wthp_pc_lst[i_chi] ** 2 * wthz_pb_lst[i_time] * np.pi ** 1.5) * \
                                  np.exp(-(vpara - ux_pb_lst[i_time]) ** 2 / wthz_pb_lst[i_time] ** 2) * \
                                  np.exp(-(vperp1 - uy_pb_lst[i_time]) ** 2 / var_wthp_pc_lst[i_chi] ** 2) * \
                                  np.exp(-(vperp2 - uz_pb_lst[i_time]) ** 2 / var_wthp_pc_lst[i_chi] ** 2)
            vdf_fit_wthz_pc_tmp = 1.e-3 * dens_pc_lst[i_time] / (
                    wthp_pc_lst[i_time] ** 2 * var_wthz_pc_lst[i_chi] * np.pi ** 1.5) * \
                                  np.exp(-(vpara - ux_pc_lst[i_time]) ** 2 / var_wthz_pc_lst[i_chi] ** 2) * \
                                  np.exp(-(vperp1 - uy_pc_lst[i_time]) ** 2 / wthp_pc_lst[i_time] ** 2) * \
                                  np.exp(-(vperp2 - uz_pc_lst[i_time]) ** 2 / wthp_pc_lst[i_time] ** 2) + \
                                  1.e-3 * dens_pb_lst[i_time] / (
                                          wthp_pb_lst[i_time] ** 2 * wthz_pb_lst[i_time] * np.pi ** 1.5) * \
                                  np.exp(-(vpara - ux_pb_lst[i_time]) ** 2 / wthz_pb_lst[i_time] ** 2) * \
                                  np.exp(-(vperp1 - uy_pb_lst[i_time]) ** 2 / wthp_pb_lst[i_time] ** 2) * \
                                  np.exp(-(vperp2 - uz_pb_lst[i_time]) ** 2 / wthp_pb_lst[i_time] ** 2)
            vdf_fit_wthz_pb_tmp = 1.e-3 * dens_pc_lst[i_time] / (
                    wthp_pc_lst[i_time] ** 2 * wthz_pc_lst[i_time] * np.pi ** 1.5) * \
                                  np.exp(-(vpara - ux_pc_lst[i_time]) ** 2 / wthz_pc_lst[i_time] ** 2) * \
                                  np.exp(-(vperp1 - uy_pc_lst[i_time]) ** 2 / wthp_pc_lst[i_time] ** 2) * \
                                  np.exp(-(vperp2 - uz_pc_lst[i_time]) ** 2 / wthp_pc_lst[i_time] ** 2) + \
                                  1.e-3 * dens_pb_lst[i_time] / (
                                          wthp_pb_lst[i_time] ** 2 * var_wthz_pb_lst[i_chi] * np.pi ** 1.5) * \
                                  np.exp(-(vpara - ux_pb_lst[i_time]) ** 2 / var_wthz_pb_lst[i_chi] ** 2) * \
                                  np.exp(-(vperp1 - uy_pb_lst[i_time]) ** 2 / wthp_pb_lst[i_time] ** 2) * \
                                  np.exp(-(vperp2 - uz_pb_lst[i_time]) ** 2 / wthp_pb_lst[i_time] ** 2)
            vdf_fit_ux_pc_tmp = 1.e-3 * dens_pc_lst[i_time] / (
                    wthp_pc_lst[i_time] ** 2 * wthz_pc_lst[i_time] * np.pi ** 1.5) * \
                                np.exp(-(vpara - var_ux_pc_lst[i_chi]) ** 2 / wthz_pc_lst[i_time] ** 2) * \
                                np.exp(-(vperp1 - uy_pc_lst[i_time]) ** 2 / wthp_pc_lst[i_time] ** 2) * \
                                np.exp(-(vperp2 - uz_pc_lst[i_time]) ** 2 / wthp_pc_lst[i_time] ** 2) + \
                                1.e-3 * dens_pb_lst[i_time] / (
                                        wthp_pb_lst[i_time] ** 2 * wthz_pb_lst[i_time] * np.pi ** 1.5) * \
                                np.exp(-(vpara - ux_pb_lst[i_time]) ** 2 / wthz_pb_lst[i_time] ** 2) * \
                                np.exp(-(vperp1 - uy_pb_lst[i_time]) ** 2 / wthp_pb_lst[i_time] ** 2) * \
                                np.exp(-(vperp2 - uz_pb_lst[i_time]) ** 2 / wthp_pb_lst[i_time] ** 2)
            vdf_fit_ux_pb_tmp = 1.e-3 * dens_pc_lst[i_time] / (
                    wthp_pc_lst[i_time] ** 2 * wthz_pc_lst[i_time] * np.pi ** 1.5) * \
                                np.exp(-(vpara - ux_pc_lst[i_time]) ** 2 / wthz_pc_lst[i_time] ** 2) * \
                                np.exp(-(vperp1 - uy_pc_lst[i_time]) ** 2 / wthp_pc_lst[i_time] ** 2) * \
                                np.exp(-(vperp2 - uz_pc_lst[i_time]) ** 2 / wthp_pc_lst[i_time] ** 2) + \
                                1.e-3 * dens_pb_lst[i_time] / (
                                        wthp_pb_lst[i_time] ** 2 * wthz_pb_lst[i_time] * np.pi ** 1.5) * \
                                np.exp(-(vpara - var_ux_pb_lst[i_chi]) ** 2 / wthz_pb_lst[i_time] ** 2) * \
                                np.exp(-(vperp1 - uy_pb_lst[i_time]) ** 2 / wthp_pb_lst[i_time] ** 2) * \
                                np.exp(-(vperp2 - uz_pb_lst[i_time]) ** 2 / wthp_pb_lst[i_time] ** 2)
            vdf_fit_uy_pc_tmp = 1.e-3 * dens_pc_lst[i_time] / (
                    wthp_pc_lst[i_time] ** 2 * wthz_pc_lst[i_time] * np.pi ** 1.5) * \
                                np.exp(-(vpara - ux_pc_lst[i_time]) ** 2 / wthz_pc_lst[i_time] ** 2) * \
                                np.exp(-(vperp1 - var_uy_pc_lst[i_chi]) ** 2 / wthp_pc_lst[i_time] ** 2) * \
                                np.exp(-(vperp2 - uz_pc_lst[i_time]) ** 2 / wthp_pc_lst[i_time] ** 2) + \
                                1.e-3 * dens_pb_lst[i_time] / (
                                        wthp_pb_lst[i_time] ** 2 * wthz_pb_lst[i_time] * np.pi ** 1.5) * \
                                np.exp(-(vpara - ux_pb_lst[i_time]) ** 2 / wthz_pb_lst[i_time] ** 2) * \
                                np.exp(-(vperp1 - uy_pb_lst[i_time]) ** 2 / wthp_pb_lst[i_time] ** 2) * \
                                np.exp(-(vperp2 - uz_pb_lst[i_time]) ** 2 / wthp_pb_lst[i_time] ** 2)
            vdf_fit_uy_pb_tmp = 1.e-3 * dens_pc_lst[i_time] / (
                    wthp_pc_lst[i_time] ** 2 * wthz_pc_lst[i_time] * np.pi ** 1.5) * \
                                np.exp(-(vpara - ux_pc_lst[i_time]) ** 2 / wthz_pc_lst[i_time] ** 2) * \
                                np.exp(-(vperp1 - uy_pc_lst[i_time]) ** 2 / wthp_pc_lst[i_time] ** 2) * \
                                np.exp(-(vperp2 - uz_pc_lst[i_time]) ** 2 / wthp_pc_lst[i_time] ** 2) + \
                                1.e-3 * dens_pb_lst[i_time] / (
                                        wthp_pb_lst[i_time] ** 2 * wthz_pb_lst[i_time] * np.pi ** 1.5) * \
                                np.exp(-(vpara - ux_pb_lst[i_time]) ** 2 / wthz_pb_lst[i_time] ** 2) * \
                                np.exp(-(vperp1 - var_uy_pb_lst[i_chi]) ** 2 / wthp_pb_lst[i_time] ** 2) * \
                                np.exp(-(vperp2 - uz_pb_lst[i_time]) ** 2 / wthp_pb_lst[i_time] ** 2)
            vdf_fit_uz_pc_tmp = 1.e-3 * dens_pc_lst[i_time] / (
                    wthp_pc_lst[i_time] ** 2 * wthz_pc_lst[i_time] * np.pi ** 1.5) * \
                                np.exp(-(vpara - ux_pc_lst[i_time]) ** 2 / wthz_pc_lst[i_time] ** 2) * \
                                np.exp(-(vperp1 - uy_pc_lst[i_time]) ** 2 / wthp_pc_lst[i_time] ** 2) * \
                                np.exp(-(vperp2 - var_uz_pc_lst[i_chi]) ** 2 / wthp_pc_lst[i_time] ** 2) + \
                                1.e-3 * dens_pb_lst[i_time] / (
                                        wthp_pb_lst[i_time] ** 2 * wthz_pb_lst[i_time] * np.pi ** 1.5) * \
                                np.exp(-(vpara - ux_pb_lst[i_time]) ** 2 / wthz_pb_lst[i_time] ** 2) * \
                                np.exp(-(vperp1 - uy_pb_lst[i_time]) ** 2 / wthp_pb_lst[i_time] ** 2) * \
                                np.exp(-(vperp2 - uz_pb_lst[i_time]) ** 2 / wthp_pb_lst[i_time] ** 2)
            vdf_fit_uz_pb_tmp = 1.e-3 * dens_pc_lst[i_time] / (
                    wthp_pc_lst[i_time] ** 2 * wthz_pc_lst[i_time] * np.pi ** 1.5) * \
                                np.exp(-(vpara - ux_pc_lst[i_time]) ** 2 / wthz_pc_lst[i_time] ** 2) * \
                                np.exp(-(vperp1 - uy_pc_lst[i_time]) ** 2 / wthp_pc_lst[i_time] ** 2) * \
                                np.exp(-(vperp2 - uz_pc_lst[i_time]) ** 2 / wthp_pc_lst[i_time] ** 2) + \
                                1.e-3 * dens_pb_lst[i_time] / (
                                        wthp_pb_lst[i_time] ** 2 * wthz_pb_lst[i_time] * np.pi ** 1.5) * \
                                np.exp(-(vpara - ux_pb_lst[i_time]) ** 2 / wthz_pb_lst[i_time] ** 2) * \
                                np.exp(-(vperp1 - uy_pb_lst[i_time]) ** 2 / wthp_pb_lst[i_time] ** 2) * \
                                np.exp(-(vperp2 - var_uz_pb_lst[i_chi]) ** 2 / wthp_pb_lst[i_time] ** 2)
            chi_n_pc_lst.append(np.nansum((vdf - vdf_fit_n_pc_tmp - mean_times_ori_sd) ** 2 / sigma))
            chi_n_pb_lst.append(np.nansum((vdf - vdf_fit_n_pb_tmp - mean_times_ori_sd) ** 2 / sigma))
            chi_wthp_pc_lst.append(np.nansum((vdf - vdf_fit_wthp_pc_tmp - mean_times_ori_sd) ** 2 / sigma))
            chi_wthp_pb_lst.append(np.nansum((vdf - vdf_fit_wthp_pb_tmp - mean_times_ori_sd) ** 2 / sigma))
            chi_wthz_pc_lst.append(np.nansum((vdf - vdf_fit_wthz_pc_tmp - mean_times_ori_sd) ** 2 / sigma))
            chi_wthz_pb_lst.append(np.nansum((vdf - vdf_fit_wthz_pb_tmp - mean_times_ori_sd) ** 2 / sigma))
            chi_ux_pc_lst.append(np.nansum((vdf - vdf_fit_ux_pc_tmp - mean_times_ori_sd) ** 2 / sigma))
            chi_ux_pb_lst.append(np.nansum((vdf - vdf_fit_ux_pb_tmp - mean_times_ori_sd) ** 2 / sigma))
            chi_uy_pc_lst.append(np.nansum((vdf - vdf_fit_uy_pc_tmp - mean_times_ori_sd) ** 2 / sigma))
            chi_uy_pb_lst.append(np.nansum((vdf - vdf_fit_uy_pb_tmp - mean_times_ori_sd) ** 2 / sigma))
            chi_uz_pc_lst.append(np.nansum((vdf - vdf_fit_uz_pc_tmp - mean_times_ori_sd) ** 2 / sigma))
            chi_uz_pb_lst.append(np.nansum((vdf - vdf_fit_uz_pb_tmp - mean_times_ori_sd) ** 2 / sigma))

    # determine the uncertainties of the parameters using (χ2 - Δχ2, χ2 + Δχ2) [confidence level: 99.$]
    # according to Numerical Recipes The Art of Scientific Computing
    delta_chi2 = 6.63
    # proton core density
    chi_fit = np.nansum((vdf - vdf_fit - mean_times_ori_sd) ** 2 / sigma)
    ind = np.argmin(np.abs(chi_n_pc_lst - chi_fit))
    min_chi, max_chi = chi_n_pc_lst[ind] - delta_chi2, chi_n_pc_lst[ind] + delta_chi2
    ind_with = np.where(np.logical_and(chi_n_pc_lst > min_chi, chi_n_pc_lst < max_chi))
    dens_pc_low_lst.append(var_n_pc_lst[np.min(ind_with)])
    dens_pc_up_lst.append(var_n_pc_lst[np.max(ind_with)])
    delta_dens_pc_low_lst.append(n_pc_chi_min - var_n_pc_lst[np.min(ind_with)])
    delta_dens_pc_up_lst.append(var_n_pc_lst[np.max(ind_with)] - n_pc_chi_min)
    # proton beam density
    chi_fit = np.nansum((vdf - vdf_fit - mean_times_ori_sd) ** 2 / sigma)
    ind = np.argmin(np.abs(chi_n_pb_lst - chi_fit))
    min_chi, max_chi = chi_n_pb_lst[ind] - delta_chi2, chi_n_pb_lst[ind] + delta_chi2
    ind_with = np.where(np.logical_and(chi_n_pb_lst > min_chi, chi_n_pb_lst < max_chi))
    dens_pb_low_lst.append(var_n_pb_lst[np.min(ind_with)])
    dens_pb_up_lst.append(var_n_pb_lst[np.max(ind_with)])
    delta_dens_pb_low_lst.append(n_pb_chi_min - var_n_pb_lst[np.min(ind_with)])
    delta_dens_pb_up_lst.append(var_n_pb_lst[np.max(ind_with)] - n_pb_chi_min)
    # proton core perpendicular thermal velocity
    chi_fit = np.nansum((vdf - vdf_fit - mean_times_ori_sd) ** 2 / sigma)
    ind = np.argmin(np.abs(chi_wthp_pc_lst - chi_fit))
    min_chi, max_chi = chi_wthp_pc_lst[ind] - delta_chi2, chi_wthp_pc_lst[ind] + delta_chi2
    ind_with = np.where(np.logical_and(chi_wthp_pc_lst > min_chi, chi_wthp_pc_lst < max_chi))
    wthp_pc_low_lst.append(var_wthp_pc_lst[np.min(ind_with)])
    wthp_pc_up_lst.append(var_wthp_pc_lst[np.max(ind_with)])
    delta_wthp_pc_low_lst.append(wthp_pc_chi_min - var_wthp_pc_lst[np.min(ind_with)])
    delta_wthp_pc_up_lst.append(var_wthp_pc_lst[np.max(ind_with)] - wthp_pc_chi_min)
    # proton beam perpendicular thermal velocity
    chi_fit = np.nansum((vdf - vdf_fit - mean_times_ori_sd) ** 2 / sigma)
    ind = np.argmin(np.abs(chi_wthp_pb_lst - chi_fit))
    min_chi, max_chi = chi_wthp_pb_lst[ind] - delta_chi2, chi_wthp_pb_lst[ind] + delta_chi2
    ind_with = np.where(np.logical_and(chi_wthp_pb_lst > min_chi, chi_wthp_pb_lst < max_chi))
    wthp_pb_low_lst.append(var_wthp_pb_lst[np.min(ind_with)])
    wthp_pb_up_lst.append(var_wthp_pb_lst[np.max(ind_with)])
    delta_wthp_pb_low_lst.append(wthp_pb_chi_min - var_wthp_pb_lst[np.min(ind_with)])
    delta_wthp_pb_up_lst.append(var_wthp_pb_lst[np.max(ind_with)] - wthp_pb_chi_min)
    # proton core parallel thermal velocity
    chi_fit = np.nansum((vdf - vdf_fit - mean_times_ori_sd) ** 2 / sigma)
    ind = np.argmin(np.abs(chi_wthz_pc_lst - chi_fit))
    min_chi, max_chi = chi_wthz_pc_lst[ind] - delta_chi2, chi_wthz_pc_lst[ind] + delta_chi2
    ind_with = np.where(np.logical_and(chi_wthz_pc_lst > min_chi, chi_wthz_pc_lst < max_chi))
    wthz_pc_low_lst.append(var_wthz_pc_lst[np.min(ind_with)])
    wthz_pc_up_lst.append(var_wthz_pc_lst[np.max(ind_with)])
    delta_wthz_pc_low_lst.append(wthz_pc_chi_min - var_wthz_pc_lst[np.min(ind_with)])
    delta_wthz_pc_up_lst.append(var_wthz_pc_lst[np.max(ind_with)] - wthz_pc_chi_min)
    # proton beam parallel thermal velocity
    chi_fit = np.nansum((vdf - vdf_fit - mean_times_ori_sd) ** 2 / sigma)
    ind = np.argmin(np.abs(chi_wthz_pb_lst - chi_fit))
    min_chi, max_chi = chi_wthz_pb_lst[ind] - delta_chi2, chi_wthz_pb_lst[ind] + delta_chi2
    ind_with = np.where(np.logical_and(chi_wthz_pb_lst > min_chi, chi_wthz_pb_lst < max_chi))
    wthz_pb_low_lst.append(var_wthz_pb_lst[np.min(ind_with)])
    wthz_pb_up_lst.append(var_wthz_pb_lst[np.max(ind_with)])
    delta_wthz_pb_low_lst.append(wthz_pb_chi_min - var_wthz_pb_lst[np.min(ind_with)])
    delta_wthz_pb_up_lst.append(var_wthz_pb_lst[np.max(ind_with)] - wthz_pb_chi_min)
    # proton core parallel velocity
    chi_fit = np.nansum((vdf - vdf_fit - mean_times_ori_sd) ** 2 / sigma)
    ind = np.argmin(np.abs(chi_ux_pc_lst - chi_fit))
    min_chi, max_chi = chi_ux_pc_lst[ind] - delta_chi2, chi_ux_pc_lst[ind] + delta_chi2
    ind_with = np.where(np.logical_and(chi_ux_pc_lst > min_chi, chi_ux_pc_lst < max_chi))
    ux_pc_low_lst.append(var_ux_pc_lst[np.min(ind_with)])
    ux_pc_up_lst.append(var_ux_pc_lst[np.max(ind_with)])
    delta_ux_pc_low_lst.append(ux_pc_chi_min - var_ux_pc_lst[np.min(ind_with)])
    delta_ux_pc_up_lst.append(var_ux_pc_lst[np.max(ind_with)] - ux_pc_chi_min)
    # proton beam parallel velocity
    chi_fit = np.nansum((vdf - vdf_fit - mean_times_ori_sd) ** 2 / sigma)
    ind = np.argmin(np.abs(chi_ux_pb_lst - chi_fit))
    min_chi, max_chi = chi_ux_pb_lst[ind] - delta_chi2, chi_ux_pb_lst[ind] + delta_chi2
    ind_with = np.where(np.logical_and(chi_ux_pb_lst > min_chi, chi_ux_pb_lst < max_chi))
    ux_pb_low_lst.append(var_ux_pb_lst[np.min(ind_with)])
    ux_pb_up_lst.append(var_ux_pb_lst[np.max(ind_with)])
    delta_ux_pb_low_lst.append(ux_pb_chi_min - var_ux_pb_lst[np.min(ind_with)])
    delta_ux_pb_up_lst.append(var_ux_pb_lst[np.max(ind_with)] - ux_pb_chi_min)
    # proton core perp1 velocity
    chi_fit = np.nansum((vdf - vdf_fit - mean_times_ori_sd) ** 2 / sigma)
    ind = np.argmin(np.abs(chi_uy_pc_lst - chi_fit))
    min_chi, max_chi = chi_uy_pc_lst[ind] - delta_chi2, chi_uy_pc_lst[ind] + delta_chi2
    ind_with = np.where(np.logical_and(chi_uy_pc_lst > min_chi, chi_uy_pc_lst < max_chi))
    uy_pc_low_lst.append(var_uy_pc_lst[np.min(ind_with)])
    uy_pc_up_lst.append(var_uy_pc_lst[np.max(ind_with)])
    delta_uy_pc_low_lst.append(uy_pc_chi_min - var_uy_pc_lst[np.min(ind_with)])
    delta_uy_pc_up_lst.append(var_uy_pc_lst[np.max(ind_with)] - uy_pc_chi_min)
    # proton beam perp1 velocity
    chi_fit = np.nansum((vdf - vdf_fit - mean_times_ori_sd) ** 2 / sigma)
    ind = np.argmin(np.abs(chi_uy_pb_lst - chi_fit))
    min_chi, max_chi = chi_uy_pb_lst[ind] - delta_chi2, chi_uy_pb_lst[ind] + delta_chi2
    ind_with = np.where(np.logical_and(chi_uy_pb_lst > min_chi, chi_uy_pb_lst < max_chi))
    uy_pb_low_lst.append(var_uy_pb_lst[np.min(ind_with)])
    uy_pb_up_lst.append(var_uy_pb_lst[np.max(ind_with)])
    delta_uy_pb_low_lst.append(uy_pb_chi_min - var_uy_pb_lst[np.min(ind_with)])
    delta_uy_pb_up_lst.append(var_uy_pb_lst[np.max(ind_with)] - uy_pb_chi_min)
    # proton core perp2 velocity
    chi_fit = np.nansum((vdf - vdf_fit - mean_times_ori_sd) ** 2 / sigma)
    ind = np.argmin(np.abs(chi_uz_pc_lst - chi_fit))
    min_chi, max_chi = chi_uz_pc_lst[ind] - delta_chi2, chi_uz_pc_lst[ind] + delta_chi2
    ind_with = np.where(np.logical_and(chi_uz_pc_lst > min_chi, chi_uz_pc_lst < max_chi))
    uz_pc_low_lst.append(var_uz_pc_lst[np.min(ind_with)])
    uz_pc_up_lst.append(var_uz_pc_lst[np.max(ind_with)])
    delta_uz_pc_low_lst.append(uz_pc_chi_min - var_uz_pc_lst[np.min(ind_with)])
    delta_uz_pc_up_lst.append(var_uz_pc_lst[np.max(ind_with)] - uz_pc_chi_min)
    # proton beam perp2 velocity
    chi_fit = np.nansum((vdf - vdf_fit - mean_times_ori_sd) ** 2 / sigma)
    ind = np.argmin(np.abs(chi_uz_pb_lst - chi_fit))
    min_chi, max_chi = chi_uz_pb_lst[ind] - delta_chi2, chi_uz_pb_lst[ind] + delta_chi2
    ind_with = np.where(np.logical_and(chi_uz_pb_lst > min_chi, chi_uz_pb_lst < max_chi))
    uz_pb_low_lst.append(var_uz_pb_lst[np.min(ind_with)])
    uz_pb_up_lst.append(var_uz_pb_lst[np.max(ind_with)])
    delta_uz_pb_low_lst.append(uz_pb_chi_min - var_uz_pb_lst[np.min(ind_with)])
    delta_uz_pb_up_lst.append(var_uz_pb_lst[np.max(ind_with)] - uz_pb_chi_min)

    # calculate errors of angle between instant magnetic field and drift velocity
    is_use_method_1_or_2_in_theta_error = 2
    theta_fit_lst.append(180 + np.rad2deg(np.arctan(
        np.sqrt((uy_pb_lst[i_time] - uy_pc_lst[i_time]) ** 2 + (uz_pb_lst[i_time] - uz_pc_lst[i_time]) ** 2) / (
                ux_pb_lst[i_time] - ux_pc_lst[i_time]))))
    if is_use_method_1_or_2_in_theta_error == 1:
        theta_B_Vd_lst_tmp = []
        for ux_pc in [ux_pc_lst[i_time], ux_pc_low_lst[i_time], ux_pc_up_lst[i_time]]:
            for uy_pc in [uy_pc_lst[i_time], uy_pc_low_lst[i_time], uy_pc_up_lst[i_time]]:
                for uz_pc in [uz_pc_lst[i_time], uz_pc_low_lst[i_time], uz_pc_up_lst[i_time]]:
                    for ux_pb in [ux_pb_lst[i_time], ux_pb_low_lst[i_time], ux_pb_up_lst[i_time]]:
                        for uy_pb in [uy_pb_lst[i_time], uy_pb_low_lst[i_time], uy_pb_up_lst[i_time]]:
                            for uz_pb in [uz_pb_lst[i_time], uz_pb_low_lst[i_time], uz_pb_up_lst[i_time]]:
                                vd_vect = [ux_pb - ux_pc, uy_pb - uy_pc, uz_pb - uz_pc]
                                theta_B_Vd_lst_tmp.append(180 +
                                                          np.rad2deg(np.arctan(
                                                              np.sqrt(vd_vect[1] ** 2 + vd_vect[2] ** 2) / vd_vect[0])))
        theta_mean_lst.append(np.mean(theta_B_Vd_lst_tmp))
        theta_std_lst.append(np.std(theta_B_Vd_lst_tmp))
        theta_ste_lst.append(np.std(theta_B_Vd_lst_tmp) / np.sqrt(len(theta_B_Vd_lst_tmp)))
    elif is_use_method_1_or_2_in_theta_error == 2:
        udx = ux_pb_lst[i_time] - ux_pc_lst[i_time]
        udy = uy_pb_lst[i_time] - uy_pc_lst[i_time]
        udz = uz_pb_lst[i_time] - uz_pc_lst[i_time]
        m = np.sqrt(udy ** 2 + udz ** 2)
        x = m / udx
        sigma_udx_low = np.sqrt(delta_ux_pc_low_lst[i_time] ** 2 + delta_ux_pb_low_lst[i_time] ** 2)
        sigma_udy_low = np.sqrt(delta_uy_pc_low_lst[i_time] ** 2 + delta_uy_pb_low_lst[i_time] ** 2)
        sigma_udz_low = np.sqrt(delta_uz_pc_low_lst[i_time] ** 2 + delta_uz_pb_low_lst[i_time] ** 2)
        sigma_m_low = np.sqrt(((udy * sigma_udy_low) ** 2 + (udz * sigma_udz_low) ** 2) / m ** 2)
        sigma_x_low = np.abs(x) * np.sqrt((sigma_m_low / m) ** 2 + (sigma_udx_low / udx) ** 2)
        theta_std_lst_low.append(np.rad2deg(sigma_x_low / (1 + x ** 2)))
        theta_ste_lst_low.append(np.rad2deg(sigma_x_low / (1 + x ** 2)))
        sigma_udx_up = np.sqrt(delta_ux_pc_up_lst[i_time] ** 2 + delta_ux_pb_up_lst[i_time] ** 2)
        sigma_udy_up = np.sqrt(delta_uy_pc_up_lst[i_time] ** 2 + delta_uy_pb_up_lst[i_time] ** 2)
        sigma_udz_up = np.sqrt(delta_uz_pc_up_lst[i_time] ** 2 + delta_uz_pb_up_lst[i_time] ** 2)
        sigma_m_up = np.sqrt(((udy * sigma_udy_up) ** 2 + (udz * sigma_udz_up) ** 2) / m ** 2)
        sigma_x_up = np.abs(x) * np.sqrt((sigma_m_up / m) ** 2 + (sigma_udx_up / udx) ** 2)
        theta_std_lst_up.append(np.rad2deg(sigma_x_up / (1 + x ** 2)))
        theta_ste_lst_up.append(np.rad2deg(sigma_x_up / (1 + x ** 2)))
        sigma_x = sigma_x_low if sigma_x_low < sigma_x_up else sigma_x_up
        theta_std_lst.append(np.rad2deg(sigma_x / (1 + x ** 2)))
        theta_ste_lst.append(np.rad2deg(sigma_x / (1 + x ** 2)))

    # plot distribution of (yi-y)/sqrt(σ) to verify the normal distribution function
    fig_dir = '/Users/psr/work/pku/nfa_beam/Figures/normal_verify/' + '-'.join([_year, _month, _day]) + '/'
    filename = 'normal_verification_' + epoch.strftime('%H%M%S.%f')[:-3] + '.png'
    fig, ax = plt.subplots(1, 1, figsize=[6, 6])
    plt.subplots_adjust(left=0.12, bottom=0.12, right=0.95, top=0.9, wspace=0.25, hspace=0.25)
    ax.plot(bin_centers, np.ma.masked_where(histogram == 0.0, histogram), marker='o', markersize=6, color='r',
            linestyle='')
    ax.plot(bin_centers, 1 / np.sqrt(2 * np.pi) * np.exp(-bin_centers ** 2 / 2), color='g', label='standard normal')
    # plt.plot(bin_centers, Gauss(bin_centers, *popt), 'r-', label='gaussian fit')
    ax.set_xlabel(r'$\frac{y_i-y_{i,fit}}{\sigma_{eff}}$', fontsize=15)
    ax.set_ylabel('PDF', fontsize=15)
    # ax.set_title(r'$\sigma^2=\sum_i\frac{(y_i-y_{i_fit})^2}{N-M}$', fontsize=15)
    ax.legend()
    # plt.show()
    plt.savefig(fig_dir + filename, dpi=300)
    plt.close()

    # plot chi-square variation
    fig_dir = '/Users/psr/work/pku/nfa_beam/Figures/chi_square/' + '-'.join([_year, _month, _day]) + '/'
    filename = 'chi_square_variation_' + epoch.strftime('%H%M%S.%f')[:-3] + '.png'
    fig, axs = plt.subplots(2, 3, figsize=[12, 8])
    plt.subplots_adjust(left=0.08, bottom=0.1, right=0.95, top=0.95, wspace=0.25, hspace=0.25)
    axs = axs.flatten()
    ax = axs[1 - 1]
    ax.plot(var_n_pc_lst, chi_n_pc_lst, label='pc', color='r')
    ax.plot(var_n_pb_lst, chi_n_pb_lst, label='pb', color='b')
    ax.axvline(x=dens_pc_lst[i_time], linestyle=':', color='r')
    ax.axvline(x=dens_pb_lst[i_time], linestyle=':', color='b')
    ax.set_xlabel(r'$N (cm^{-3})$', fontsize=15)
    ax.set_ylabel(r'$\chi^2$', fontsize=15)
    ax.legend()
    ax = axs[2 - 1]
    ax.plot(var_wthp_pc_lst, chi_wthp_pc_lst, color='r')
    ax.plot(var_wthp_pb_lst, chi_wthp_pb_lst, color='b')
    ax.axvline(x=wthp_pc_lst[i_time], linestyle=':', color='r')
    ax.axvline(x=wthp_pb_lst[i_time], linestyle=':', color='b')
    ax.set_xlabel(r'$w_{th,\perp} (km/s)$', fontsize=15)
    ax.set_ylabel(r'$\chi^2$', fontsize=15)
    ax = axs[3 - 1]
    ax.plot(var_wthz_pc_lst, chi_wthz_pc_lst, color='r')
    ax.plot(var_wthz_pb_lst, chi_wthz_pb_lst, color='b')
    ax.axvline(x=wthz_pc_lst[i_time], linestyle=':', color='r')
    ax.axvline(x=wthz_pb_lst[i_time], linestyle=':', color='b')
    ax.set_xlabel(r'$w_{th,||} (km/s)$', fontsize=15)
    ax.set_ylabel(r'$\chi^2$', fontsize=15)
    ax = axs[4 - 1]
    ax.plot(var_ux_pc_lst, chi_ux_pc_lst, color='r')
    ax.plot(var_ux_pb_lst, chi_ux_pb_lst, color='b')
    ax.axvline(x=ux_pc_lst[i_time], linestyle=':', color='r')
    ax.axvline(x=ux_pb_lst[i_time], linestyle=':', color='b')
    ax.set_xlabel(r'$U_{||} (km/s)$', fontsize=15)
    ax.set_ylabel(r'$\chi^2$', fontsize=15)
    ax = axs[5 - 1]
    ax.plot(var_uy_pc_lst, chi_uy_pc_lst, color='r')
    ax.plot(var_uy_pb_lst, chi_uy_pb_lst, color='b')
    ax.axvline(x=uy_pc_lst[i_time], linestyle=':', color='r')
    ax.axvline(x=uy_pb_lst[i_time], linestyle=':', color='b')
    ax.set_xlabel(r'$U_{\perp,1} (km/s)$', fontsize=15)
    ax.set_ylabel(r'$\chi^2$', fontsize=15)
    ax = axs[6 - 1]
    ax.plot(var_uz_pc_lst, chi_uz_pc_lst, color='r')
    ax.plot(var_uz_pb_lst, chi_uz_pb_lst, color='b')
    ax.axvline(x=uz_pc_lst[i_time], linestyle=':', color='r')
    ax.axvline(x=uz_pb_lst[i_time], linestyle=':', color='b')
    ax.set_xlabel(r'$U_{\perp,2} (km/s)$', fontsize=15)
    ax.set_ylabel(r'$\chi^2$', fontsize=15)
    # plt.show()
    plt.savefig(fig_dir + filename, dpi=300)
    plt.close()

    fig_dir = '/Users/psr/work/pku/nfa_beam/Figures/chi_square/' + '-'.join([_year, _month, _day]) + '/'
    filename = 'chi_square_variation_rtn_bulk_velocity_' + epoch.strftime('%H%M%S.%f')[:-3] + '.png'
    fig, axs = plt.subplots(2, 3, figsize=[12, 8])
    plt.subplots_adjust(left=0.08, bottom=0.1, right=0.95, top=0.95, wspace=0.25, hspace=0.25)
    axs = axs.flatten()
    ax = axs[1 - 1]
    ax.plot(var_ur_pc_lst, chi_ur_pc_lst, label='pc', color='r')
    ax.axvline(x=ur_pc_lst[i_time], linestyle=':', color='r')
    ax.set_xlabel(r'$U_{pc,R} (km/s)$', fontsize=15)
    ax.set_ylabel(r'$\chi^2$', fontsize=15)
    ax = axs[2 - 1]
    ax.plot(var_ut_pc_lst, chi_ut_pc_lst, label='pc', color='r')
    ax.axvline(x=ut_pc_lst[i_time], linestyle=':', color='r')
    ax.set_xlabel(r'$U_{pc,T} (km/s)$', fontsize=15)
    ax.set_ylabel(r'$\chi^2$', fontsize=15)
    ax = axs[3 - 1]
    ax.plot(var_un_pc_lst, chi_un_pc_lst, label='pc', color='r')
    ax.axvline(x=un_pc_lst[i_time], linestyle=':', color='r')
    ax.set_xlabel(r'$U_{pc,N} (km/s)$', fontsize=15)
    ax.set_ylabel(r'$\chi^2$', fontsize=15)
    ax = axs[4 - 1]
    ax.plot(var_ur_pb_lst, chi_ur_pb_lst, label='pc', color='b')
    ax.axvline(x=ur_pb_lst[i_time], linestyle=':', color='b')
    ax.set_xlabel(r'$U_{pb,R} (km/s)$', fontsize=15)
    ax.set_ylabel(r'$\chi^2$', fontsize=15)
    ax = axs[5 - 1]
    ax.plot(var_ut_pb_lst, chi_ut_pb_lst, label='pc', color='b')
    ax.axvline(x=ut_pb_lst[i_time], linestyle=':', color='b')
    ax.set_xlabel(r'$U_{pb,T} (km/s)$', fontsize=15)
    ax.set_ylabel(r'$\chi^2$', fontsize=15)
    ax = axs[6 - 1]
    ax.plot(var_un_pb_lst, chi_un_pb_lst, label='pc', color='b')
    ax.axvline(x=un_pb_lst[i_time], linestyle=':', color='b')
    ax.set_xlabel(r'$U_{pb,N} (km/s)$', fontsize=15)
    ax.set_ylabel(r'$\chi^2$', fontsize=15)
    # plt.show()
    plt.savefig(fig_dir + filename, dpi=300)
    plt.close()
plot_confidence_level()

time_beg = str(np.min(hour_lst)) + str(np.min(min_lst)) + str(np.min(sec_lst))
time_end = str(np.max(hour_lst)) + str(np.max(min_lst)) + str(np.max(sec_lst))
if is_use_method_1_or_2_in_theta_error == 2:
    print('theta standard deviations:')
    print('fit value (average): ', np.mean(theta_fit_lst))
    print('mean value (low, up): ', np.mean(theta_std_lst_low), ',', np.mean(theta_std_lst_up))
    print('median value (low, up): ', np.median(theta_std_lst_low), ',', np.median(theta_std_lst_up))
    dir_solution = '/Users/psr/work/pku/nfa_beam/Data/theta/'
    filename = 'theta_fit_sigma_low_up' + time_beg + '_' + time_end + '.txt'
    solution = np.hstack([np.array(theta_fit_lst).reshape(-1, 1), np.array(theta_std_lst_low).reshape(-1, 1),
                          np.array(theta_std_lst_up).reshape(-1, 1)])
    np.savetxt(dir_solution + filename, solution, fmt='%.10f')

# plot time series of the fitting parameters
fig_dir = '/Users/psr/work/pku/nfa_beam/Figures/time_series/' + '-'.join([_year, _month, _day]) + '/'
filename = 'polyfill_time_series_uncertainty_' + time_beg + '_' + time_end + '.png'
fig, axs = plt.subplots(6, 1, figsize=[10, 10], sharex=True)
plt.subplots_adjust(left=0.12, bottom=0.1, right=0.88, top=0.95, wspace=0.2, hspace=0.1)
ax = axs[0]
y_low = np.array(ux_pc_low_lst)
y_up = np.array(ux_pc_up_lst)
ax.fill_between(datetime_lst, y_low, y_up, color='palegreen')
ax.plot(datetime_lst, ux_pc_lst, color='k')
ax.set_ylabel(r'$U_{pc,||} (km/s)$', fontsize=10)

ax = axs[1]
y_low = np.array(ux_pb_low_lst)
y_up = np.array(ux_pb_up_lst)
ax.fill_between(datetime_lst, y_low, y_up, color='palegreen')
ax.plot(datetime_lst, ux_pb_lst, color='k')
ax.set_ylabel(r'$U_{pb,||} (km/s)$', fontsize=10)

ax = axs[2]
y_low = np.array(uy_pc_low_lst)
y_up = np.array(uy_pc_up_lst)
ax.fill_between(datetime_lst, y_low, y_up, color='palegreen')
ax.plot(datetime_lst, uy_pc_lst, color='k')
ax.set_ylabel(r'$U_{pc,\perp 1} (km/s)$', fontsize=10)

ax = axs[3]
y_low = np.array(uy_pb_low_lst)
y_up = np.array(uy_pb_up_lst)
ax.fill_between(datetime_lst, y_low, y_up, color='palegreen')
ax.plot(datetime_lst, uy_pb_lst, color='k')
ax.set_ylabel(r'$U_{pb,\perp 1} (km/s)$', fontsize=10)

ax = axs[4]
y_low = np.array(uz_pc_low_lst)
y_up = np.array(uz_pc_up_lst)
ax.fill_between(datetime_lst, y_low, y_up, color='palegreen')
ax.plot(datetime_lst, uz_pc_lst, color='k')
ax.set_ylabel(r'$U_{pc,\perp 2} (km/s)$', fontsize=10)

ax = axs[5]
y_low = np.array(uz_pb_low_lst)
y_up = np.array(uz_pb_up_lst)
ax.fill_between(datetime_lst, y_low, y_up, color='palegreen')
ax.plot(datetime_lst, uz_pb_lst, color='k')
ax.set_ylabel(r'$U_{pb,\perp 1} (km/s)$', fontsize=10)
plt.savefig(fig_dir + filename, dpi=300)

filename = 'polyfill_time_series_uncertainty_rtn_bulk_velocity_' + time_beg + '_' + time_end + '.png'
fig, axs = plt.subplots(7, 1, figsize=[10, 10], sharex=True)
plt.subplots_adjust(left=0.12, bottom=0.1, right=0.88, top=0.95, wspace=0.2, hspace=0.1)
ax = axs[0]
# y_low = np.array(ur_pc_low_lst).reshape(-1)
# y_up = np.array(ur_pc_up_lst).reshape(-1)
y_low = np.array(ur_pc_lst).reshape(-1) - np.array(delta_ur_pc_low_lst).reshape(-1)
y_up = np.array(ur_pc_lst).reshape(-1) + np.array(delta_ur_pc_up_lst).reshape(-1)
ax.fill_between(datetime_lst, y_low, y_up, color='palegreen')
ax.plot(datetime_lst, ur_pc_lst, color='k')
ax.set_ylabel(r'$U_{pc,R} (km/s)$', fontsize=12)

ax = axs[1]
# y_low = np.array(ur_pb_low_lst).reshape(-1)
# y_up = np.array(ur_pb_up_lst).reshape(-1)
y_low = np.array(ur_pb_lst).reshape(-1) - np.array(delta_ur_pb_low_lst).reshape(-1)
y_up = np.array(ur_pb_lst).reshape(-1) + np.array(delta_ur_pb_up_lst).reshape(-1)
ax.fill_between(datetime_lst, y_low, y_up, color='palegreen')
ax.plot(datetime_lst, ur_pb_lst, color='k')
ax.set_ylabel(r'$U_{pb,R} (km/s)$', fontsize=12)

ax = axs[2]
# y_low = np.array(ut_pc_low_lst).reshape(-1)
# y_up = np.array(ut_pc_up_lst).reshape(-1)
y_low = np.array(ut_pc_lst).reshape(-1) - np.array(delta_ut_pc_low_lst).reshape(-1)
y_up = np.array(ut_pc_lst).reshape(-1) + np.array(delta_ut_pc_up_lst).reshape(-1)
ax.fill_between(datetime_lst, y_low, y_up, color='palegreen')
ax.plot(datetime_lst, ut_pc_lst, color='k')
ax.set_ylabel(r'$U_{pc,T} (km/s)$', fontsize=12)

ax = axs[3]
# y_low = np.array(ut_pb_low_lst).reshape(-1)
# y_up = np.array(ut_pb_up_lst).reshape(-1)
y_low = np.array(ut_pb_lst).reshape(-1) - np.array(delta_ut_pb_low_lst).reshape(-1)
y_up = np.array(ut_pb_lst).reshape(-1) + np.array(delta_ut_pb_up_lst).reshape(-1)
ax.fill_between(datetime_lst, y_low, y_up, color='palegreen')
ax.plot(datetime_lst, ut_pb_lst, color='k')
ax.set_ylabel(r'$U_{pb,T} (km/s)$', fontsize=12)

ax = axs[4]
# y_low = np.array(un_pc_low_lst).reshape(-1)
# y_up = np.array(un_pc_up_lst).reshape(-1)
y_low = np.array(un_pc_lst).reshape(-1) - np.array(delta_un_pc_low_lst).reshape(-1)
y_up = np.array(un_pc_lst).reshape(-1) + np.array(delta_un_pc_up_lst).reshape(-1)
ax.fill_between(datetime_lst, y_low, y_up, color='palegreen')
ax.plot(datetime_lst, un_pc_lst, color='k')
ax.set_ylabel(r'$U_{pc,N} (km/s)$', fontsize=12)

ax = axs[5]
# y_low = np.array(un_pb_low_lst).reshape(-1)
# y_up = np.array(un_pb_up_lst).reshape(-1)
y_low = np.array(un_pb_lst).reshape(-1) - np.array(delta_un_pb_low_lst).reshape(-1)
y_up = np.array(un_pb_lst).reshape(-1) + np.array(delta_un_pb_up_lst).reshape(-1)
ax.fill_between(datetime_lst, y_low, y_up, color='palegreen')
ax.plot(datetime_lst, un_pb_lst, color='k')
ax.set_ylabel(r'$U_{pb,N} (km/s)$', fontsize=12)

ax = axs[6]
y_low = np.array(theta_fit_lst).reshape(-1) - np.array(theta_ste_lst).reshape(-1)
y_up = np.array(theta_fit_lst).reshape(-1) + np.array(theta_ste_lst).reshape(-1)
ax.fill_between(datetime_lst, y_low, y_up, color='palegreen')
ax.plot(datetime_lst, theta_fit_lst, color='k')
ax.set_ylabel(r'$\theta_{B_0,V_d} (\circ) $', fontsize=12)
ax.set_ylim([90, 180])
ax.yaxis.set_major_locator(MultipleLocator(30))
ax.xaxis.set_major_locator(mdates.MinuteLocator(byminute=range(0, 60, 5)))
ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
plt.savefig(fig_dir + filename, dpi=300)

filename = 'time_series_' + time_beg + '_' + time_end + '.png'
fig, axs = plt.subplots(6, 1, figsize=[10, 8], sharex=True)
plt.subplots_adjust(left=0.12, bottom=0.1, right=0.88, top=0.95, wspace=0.2, hspace=0.1)
ax = axs[1 - 1]
ax.plot_date(datetime_lst, dens_pc_lst, 'r', label='pc')
ax.plot_date(datetime_lst, dens_pb_lst, 'b', label='pb')
ax.set_ylabel(r'$n (cm^{-3})$', fontsize=10)
ax.legend()
ax = axs[2 - 1]
ax.plot_date(datetime_lst, wthp_pc_lst, 'r', label=r'$w_{th,\perp}$')
ax.plot_date(datetime_lst, wthz_pc_lst, 'b', label=r'$w_{th,||}$')
ax.set_ylabel(r'$w_{th,pc} (km/s)$', fontsize=10)
ax.legend()
ax = axs[3 - 1]
ax.plot_date(datetime_lst, wthp_pb_lst, 'r', label=r'$w_{th,\perp}$')
ax.plot_date(datetime_lst, wthz_pb_lst, 'b', label=r'$w_{th,||}$')
ax.set_ylabel(r'$w_{th,pb} (km/s)$', fontsize=10)
ax.legend()
ax = axs[4 - 1]
ax.plot_date(datetime_lst, ux_pc_lst, 'r', label='pc')
ax.plot_date(datetime_lst, ux_pb_lst, 'b', label='pb')
ax.set_ylabel(r'$U_x (km/s)$', fontsize=10)
ax.legend()
ax = axs[5 - 1]
ax.plot_date(datetime_lst, uy_pc_lst, 'r', label='pc')
ax.plot_date(datetime_lst, uy_pb_lst, 'b', label='pb')
ax.set_ylabel(r'$U_y (km/s)$', fontsize=10)
ax.legend()
ax = axs[6 - 1]
ax.plot_date(datetime_lst, uz_pc_lst, 'r', label='pc')
ax.plot_date(datetime_lst, uz_pb_lst, 'b', label='pb')
ax.set_ylabel(r'$U_z (km/s)$', fontsize=10)
ax.legend()
plt.savefig(fig_dir + filename, dpi=300)

# plot time series of the fitting parameters (with uncertainty)
fig_dir = '/Users/psr/work/pku/nfa_beam/Figures/time_series/' + '-'.join([_year, _month, _day]) + '/'
filename = 'time_series_uncertainty' + time_beg + '_' + time_end + '.png'
fig, axs = plt.subplots(6, 2, figsize=[16, 10], sharex=True)
plt.subplots_adjust(left=0.12, bottom=0.1, right=0.88, top=0.95, wspace=0.2, hspace=0.1)
ax = axs[0, 0]
yerr_low_68pct = np.array(dens_pc_lst) - np.array(dens_pc_low_lst_68pct)
yerr_up_68pct = np.array(dens_pc_up_lst_68pct) - np.array(dens_pc_lst)
yerr_low_90pct = np.array(dens_pc_lst) - np.array(dens_pc_low_lst_90pct)
yerr_up_90pct = np.array(dens_pc_up_lst_90pct) - np.array(dens_pc_lst)
err_arr1 = np.vstack([yerr_low_90pct, yerr_up_90pct])
err_arr2 = np.vstack([yerr_low_68pct, yerr_up_68pct])
ax.errorbar(datetime_lst, dens_pc_lst, yerr=err_arr1, color='g', ls='none')
ax.errorbar(datetime_lst, dens_pc_lst, yerr=err_arr2, color='b', ls='none')
ax.plot_date(datetime_lst, dens_pc_lst, 'k', marker='o', markersize=2)
ax.set_ylabel(r'$n_{pc} (cm^{-3})$', fontsize=10)

ax = axs[1, 0]
yerr_low_68pct = np.array(dens_pb_lst) - np.array(dens_pb_low_lst_68pct)
yerr_up_68pct = np.array(dens_pb_up_lst_68pct) - np.array(dens_pb_lst)
yerr_low_90pct = np.array(dens_pb_lst) - np.array(dens_pb_low_lst_90pct)
yerr_up_90pct = np.array(dens_pb_up_lst_90pct) - np.array(dens_pb_lst)
err_arr1 = np.vstack([yerr_low_90pct, yerr_up_90pct])
err_arr2 = np.vstack([yerr_low_68pct, yerr_up_68pct])
ax.errorbar(datetime_lst, dens_pb_lst, yerr=err_arr1, color='g', ls='none')
ax.errorbar(datetime_lst, dens_pb_lst, yerr=err_arr2, color='b', ls='none')
ax.plot_date(datetime_lst, dens_pb_lst, 'k', marker='o', markersize=2)
ax.set_ylabel(r'$n_{pb} (cm^{-3})$', fontsize=10)

ax = axs[2, 0]
yerr_low_68pct = np.array(wthp_pc_lst) - np.array(wthp_pc_low_lst_68pct)
yerr_up_68pct = np.array(wthp_pc_up_lst_68pct) - np.array(wthp_pc_lst)
yerr_low_90pct = np.array(wthp_pc_lst) - np.array(wthp_pc_low_lst_90pct)
yerr_up_90pct = np.array(wthp_pc_up_lst_90pct) - np.array(wthp_pc_lst)
err_arr1 = np.vstack([yerr_low_90pct, yerr_up_90pct])
err_arr2 = np.vstack([yerr_low_68pct, yerr_up_68pct])
ax.errorbar(datetime_lst, wthp_pc_lst, yerr=err_arr1, color='g', ls='none')
ax.errorbar(datetime_lst, wthp_pc_lst, yerr=err_arr2, color='b', ls='none')
ax.plot_date(datetime_lst, wthp_pc_lst, 'k', marker='o', markersize=2)
ax.set_ylabel(r'$w_{thp,pc} (km/s)$', fontsize=10)

ax = axs[3, 0]
yerr_low_68pct = np.array(wthp_pb_lst) - np.array(wthp_pb_low_lst_68pct)
yerr_up_68pct = np.array(wthp_pb_up_lst_68pct) - np.array(wthp_pb_lst)
yerr_low_90pct = np.array(wthp_pb_lst) - np.array(wthp_pb_low_lst_90pct)
yerr_up_90pct = np.array(wthp_pb_up_lst_90pct) - np.array(wthp_pb_lst)
err_arr1 = np.vstack([yerr_low_90pct, yerr_up_90pct])
err_arr2 = np.vstack([yerr_low_68pct, yerr_up_68pct])
ax.errorbar(datetime_lst, wthp_pb_lst, yerr=err_arr1, color='g', ls='none')
ax.errorbar(datetime_lst, wthp_pb_lst, yerr=err_arr2, color='b', ls='none')
ax.plot_date(datetime_lst, wthp_pb_lst, 'k', marker='o', markersize=2)
ax.set_ylabel(r'$w_{thp,pb} (km/s)$', fontsize=10)

ax = axs[4, 0]
yerr_low_68pct = np.array(wthz_pc_lst) - np.array(wthz_pc_low_lst_68pct)
yerr_up_68pct = np.array(wthz_pc_up_lst_68pct) - np.array(wthz_pc_lst)
yerr_low_90pct = np.array(wthz_pc_lst) - np.array(wthz_pc_low_lst_90pct)
yerr_up_90pct = np.array(wthz_pc_up_lst_90pct) - np.array(wthz_pc_lst)
err_arr1 = np.vstack([yerr_low_90pct, yerr_up_90pct])
err_arr2 = np.vstack([yerr_low_68pct, yerr_up_68pct])
ax.errorbar(datetime_lst, wthz_pc_lst, yerr=err_arr1, color='g', ls='none')
ax.errorbar(datetime_lst, wthz_pc_lst, yerr=err_arr2, color='b', ls='none')
ax.plot_date(datetime_lst, wthz_pc_lst, 'k', marker='o', markersize=2)
ax.set_ylabel(r'$w_{thz,pc} (km/s)$', fontsize=10)

ax = axs[5, 0]
yerr_low_68pct = np.array(wthz_pb_lst) - np.array(wthz_pb_low_lst_68pct)
yerr_up_68pct = np.array(wthz_pb_up_lst_68pct) - np.array(wthz_pb_lst)
yerr_low_90pct = np.array(wthz_pb_lst) - np.array(wthz_pb_low_lst_90pct)
yerr_up_90pct = np.array(wthz_pb_up_lst_90pct) - np.array(wthz_pb_lst)
err_arr1 = np.vstack([yerr_low_90pct, yerr_up_90pct])
err_arr2 = np.vstack([yerr_low_68pct, yerr_up_68pct])
ax.errorbar(datetime_lst, wthz_pb_lst, yerr=err_arr1, color='g', ls='none')
ax.errorbar(datetime_lst, wthz_pb_lst, yerr=err_arr2, color='b', ls='none')
ax.plot_date(datetime_lst, wthz_pb_lst, 'k', marker='o', markersize=2)
ax.set_ylabel(r'$w_{thz,pb} (km/s)$', fontsize=10)

ax = axs[0, 1]
yerr_low_68pct = np.array(ux_pc_lst) - np.array(ux_pc_low_lst_68pct)
yerr_up_68pct = np.array(ux_pc_up_lst_68pct) - np.array(ux_pc_lst)
yerr_low_90pct = np.array(ux_pc_lst) - np.array(ux_pc_low_lst_90pct)
yerr_up_90pct = np.array(ux_pc_up_lst_90pct) - np.array(ux_pc_lst)
err_arr1 = np.vstack([yerr_low_90pct, yerr_up_90pct])
err_arr2 = np.vstack([yerr_low_68pct, yerr_up_68pct])
ax.errorbar(datetime_lst, ux_pc_lst, yerr=err_arr1, color='g', ls='none')
ax.errorbar(datetime_lst, ux_pc_lst, yerr=err_arr2, color='b', ls='none')
ax.plot_date(datetime_lst, ux_pc_lst, 'k', marker='o', markersize=2)
ax.set_ylabel(r'$U_{x,pc} (km/s)$', fontsize=10)

ax = axs[1, 1]
yerr_low_68pct = np.array(ux_pb_lst) - np.array(ux_pb_low_lst_68pct)
yerr_up_68pct = np.array(ux_pb_up_lst_68pct) - np.array(ux_pb_lst)
yerr_low_90pct = np.array(ux_pb_lst) - np.array(ux_pb_low_lst_90pct)
yerr_up_90pct = np.array(ux_pb_up_lst_90pct) - np.array(ux_pb_lst)
err_arr1 = np.vstack([yerr_low_90pct, yerr_up_90pct])
err_arr2 = np.vstack([yerr_low_68pct, yerr_up_68pct])
ax.errorbar(datetime_lst, ux_pb_lst, yerr=err_arr1, color='g', ls='none')
ax.errorbar(datetime_lst, ux_pb_lst, yerr=err_arr2, color='b', ls='none')
ax.plot_date(datetime_lst, ux_pb_lst, 'k', marker='o', markersize=2)
ax.set_ylabel(r'$U_{x,pb} (km/s)$', fontsize=10)

ax = axs[2, 1]
yerr_low_68pct = np.array(uy_pc_lst) - np.array(uy_pc_low_lst_68pct)
yerr_up_68pct = np.array(uy_pc_up_lst_68pct) - np.array(uy_pc_lst)
yerr_low_90pct = np.array(uy_pc_lst) - np.array(uy_pc_low_lst_90pct)
yerr_up_90pct = np.array(uy_pc_up_lst_90pct) - np.array(uy_pc_lst)
err_arr1 = np.vstack([yerr_low_90pct, yerr_up_90pct])
err_arr2 = np.vstack([yerr_low_68pct, yerr_up_68pct])
ax.errorbar(datetime_lst, uy_pc_lst, yerr=err_arr1, color='g', ls='none')
ax.errorbar(datetime_lst, uy_pc_lst, yerr=err_arr2, color='b', ls='none')
ax.plot_date(datetime_lst, uy_pc_lst, 'k', marker='o', markersize=2)
ax.set_ylabel(r'$U_{x,pc} (km/s)$', fontsize=10)

ax = axs[3, 1]
yerr_low_68pct = np.array(uy_pb_lst) - np.array(uy_pb_low_lst_68pct)
yerr_up_68pct = np.array(uy_pb_up_lst_68pct) - np.array(uy_pb_lst)
yerr_low_90pct = np.array(uy_pb_lst) - np.array(uy_pb_low_lst_90pct)
yerr_up_90pct = np.array(uy_pb_up_lst_90pct) - np.array(uy_pb_lst)
err_arr1 = np.vstack([yerr_low_90pct, yerr_up_90pct])
err_arr2 = np.vstack([yerr_low_68pct, yerr_up_68pct])
ax.errorbar(datetime_lst, uy_pb_lst, yerr=err_arr1, color='g', ls='none')
ax.errorbar(datetime_lst, uy_pb_lst, yerr=err_arr2, color='b', ls='none')
ax.plot_date(datetime_lst, uy_pb_lst, 'k', marker='o', markersize=2)
ax.set_ylabel(r'$U_{x,pb} (km/s)$', fontsize=10)

ax = axs[4, 1]
yerr_low_68pct = np.array(ux_pc_lst) - np.array(ux_pc_low_lst_68pct)
yerr_up_68pct = np.array(ux_pc_up_lst_68pct) - np.array(ux_pc_lst)
yerr_low_90pct = np.array(uz_pc_lst) - np.array(uz_pc_low_lst_90pct)
yerr_up_90pct = np.array(uz_pc_up_lst_90pct) - np.array(uz_pc_lst)
err_arr1 = np.vstack([yerr_low_90pct, yerr_up_90pct])
err_arr2 = np.vstack([yerr_low_68pct, yerr_up_68pct])
ax.errorbar(datetime_lst, uz_pc_lst, yerr=err_arr1, color='g', ls='none')
ax.errorbar(datetime_lst, uz_pc_lst, yerr=err_arr2, color='b', ls='none')
ax.plot_date(datetime_lst, uz_pc_lst, 'k', marker='o', markersize=2)
ax.set_ylabel(r'$U_{z,pc} (km/s)$', fontsize=10)

ax = axs[5, 1]
yerr_low_68pct = np.array(uz_pb_lst) - np.array(uz_pb_low_lst_68pct)
yerr_up_68pct = np.array(uz_pb_up_lst_68pct) - np.array(uz_pb_lst)
yerr_low_90pct = np.array(uz_pb_lst) - np.array(uz_pb_low_lst_90pct)
yerr_up_90pct = np.array(uz_pb_up_lst_90pct) - np.array(uz_pb_lst)
err_arr1 = np.vstack([yerr_low_90pct, yerr_up_90pct])
err_arr2 = np.vstack([yerr_low_68pct, yerr_up_68pct])
ax.errorbar(datetime_lst, uz_pb_lst, yerr=err_arr1, color='g', ls='none')
ax.errorbar(datetime_lst, uz_pb_lst, yerr=err_arr2, color='b', ls='none')
ax.plot_date(datetime_lst, uz_pb_lst, 'k', marker='o', markersize=2)
ax.set_ylabel(r'$U_{z,pb} (km/s)$', fontsize=10)
plt.savefig(fig_dir + filename, dpi=300)


# filename = 'time_series_bulk_veclocity_rtn_' + time_beg + '_' + time_end + '.png'
fig, axs = plt.subplots(6, 1, figsize=[10, 8], sharex=True)
plt.subplots_adjust(left=0.12, bottom=0.1, right=0.88, top=0.95, wspace=0.2, hspace=0.1)
ax = axs[1 - 1]
ax.plot_date(datetime_lst, ur_pc_lst, 'r', label='pc')
ax.set_ylabel(r'$U_{pc,R} (km/s)$', fontsize=10)
ax = axs[2 - 1]
ax.plot_date(datetime_lst, ur_pb_lst, 'b', label='pb')
ax.set_ylabel(r'$U_{pb,R} (km/s)$', fontsize=10)
ax = axs[3 - 1]
ax.plot_date(datetime_lst, ut_pc_lst, 'r', label='pc')
ax.set_ylabel(r'$U_{pc,T} (km/s)$', fontsize=10)
ax = axs[4 - 1]
ax.plot_date(datetime_lst, ut_pb_lst, 'b', label='pb')
ax.set_ylabel(r'$U_{pb,T} (km/s)$', fontsize=10)
ax = axs[5 - 1]
ax.plot_date(datetime_lst, un_pc_lst, 'r', label='pc')
ax.set_ylabel(r'$U_{pc,N} (km/s)$', fontsize=10)
ax = axs[6 - 1]
ax.plot_date(datetime_lst, un_pb_lst, 'b', label='pb')
ax.set_ylabel(r'$U_{pb,N} (km/s)$', fontsize=10)
plt.savefig(fig_dir + filename, dpi=300)
