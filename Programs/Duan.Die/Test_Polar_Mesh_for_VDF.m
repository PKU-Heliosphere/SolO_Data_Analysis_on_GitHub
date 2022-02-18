%% This script is used to plot SolO PAS high resolution and fitting in RTN coordinates
%close all;
%clc; clear all; %#ok<CLALL>
event_type = 'magnetic reconnection'
onesec = datenum(2010,1,1,1,1,1)-datenum(2010,1,1,1,1,0);
m_p = 1.673e-27;%kg
eV = 1.602e-19;%J
date_str = '20201014';
hour_beg = 22; min_beg = 45; sec_beg = 00;
hour_end = 23; min_end = 05; sec_end = 00;
year_str = date_str(1:4); mon_str = date_str(5:6); day_str = date_str(7:8);
year = str2num(year_str); month = str2num(mon_str); day = str2num(day_str);
time_beg = datenum(year,month,day,hour_beg,min_beg,sec_beg);
time_end = datenum(year,month,day,hour_end,min_end,sec_end);
mag_hour_beg = 22; mag_min_beg = 45; mag_sec_beg = 00;
mag_hour_end = 23; mag_min_end = 05; mag_sec_end = 00;
mag_time_beg = datenum(year,month,day,mag_hour_beg,mag_min_beg,mag_sec_beg);
mag_time_end = datenum(year,month,day,mag_hour_end,mag_min_end,mag_sec_end);
savedir = 'Figures\Reconnection\'
%% Load mag data
magdir =  ['D:\SolOData\'];
magdir = [magdir 'solo_l2_mag-rtn-normal_' date_str '_v01.cdf'];
epochmag = spdfcdfread(magdir,'Variables',{'EPOCH'});
B_rtn = double(spdfcdfread(magdir,'Variables',{'B_RTN'}));
sub_magplot = find(epochmag >= mag_time_beg & epochmag <= mag_time_end);
epochmag_plot = epochmag(sub_magplot);
br_plot = B_rtn(sub_magplot,1);
bt_plot = B_rtn(sub_magplot,2);
bn_plot = B_rtn(sub_magplot,3);
%% Calc and Plot MVA
[BL,BM,BN] = calc_mag_MVA(br_plot,bt_plot,bn_plot);
BLMN = [BL;BM;BN]*[br_plot bt_plot bn_plot]';
figure;
subplot(2,1,1)
plot(epochmag_plot,[br_plot bt_plot bn_plot])
hold on
plot(epochmag_plot,vecnorm([br_plot bt_plot bn_plot],2,2))
datetick('x','HH:MM')
legend('B_R','B_T','B_N','|B|')
subplot(2,1,2)
plot(epochmag_plot,BLMN)
hold on
plot(epochmag_plot,vecnorm([br_plot bt_plot bn_plot],2,2))
datetick('x','HH:MM')
legend('B_L','B_M','B_N','|B|')
%% load moments
momsdir = ['D:\SolOData\solo_l2_swa-pas-grnd-mom_' date_str '_v01.cdf'];
epochmoms = spdfcdfread(momsdir,'Variables',{'Epoch'});
N_rtn = double(spdfcdfread(momsdir,'Variables',{'N'}));
T_rtn = double(spdfcdfread(momsdir,'Variables',{'T'}));
V_rtn = double(spdfcdfread(momsdir,'Variables',{'V_RTN'}));
%%
epochmoms_plot = epochmoms(epochmoms>mag_time_beg & epochmoms<mag_time_end);
N_rtn_plot = N_rtn(epochmoms>mag_time_beg & epochmoms<mag_time_end);
T_rtn_plot = T_rtn(epochmoms>mag_time_beg & epochmoms<mag_time_end);
V_rtn_plot = V_rtn(epochmoms>mag_time_beg & epochmoms<mag_time_end,:);
VLMN_moms = [BL;BM;BN]*V_rtn_plot';
BLMN_moms = interp1(epochmag_plot,BLMN',epochmoms_plot);
v_A_arr = interp1(epochmag_plot,vecnorm([br_plot bt_plot bn_plot],2,2),epochmoms_plot)*1e-9./sqrt(m_p.*N_rtn_plot*1e6*4e-7*pi);
Omega_p = eV*interp1(epochmag_plot,vecnorm([br_plot bt_plot bn_plot],2,2),epochmoms_plot)*1e-9/m_p;
di_arr = v_A_arr./Omega_p;
beta_arr = N_rtn_plot*1e6.*T_rtn_plot*eV./(vecnorm(BLMN_moms,2,2)*1e-9).^2*2*4e-7*pi;

%% load vdf
cdfdir = ['D:\SolOData\solo_l2_swa-pas-vdf_' date_str '_v01.cdf'];
epoch = spdfcdfread(cdfdir,'Variables',{'Epoch'});
energybin = double(spdfcdfread(cdfdir,'Variables',{'Energy'}));
azibin = double(spdfcdfread(cdfdir,'Variables',{'Azimuth'}));
elebin = double(spdfcdfread(cdfdir,'Variables',{'Elevation'}));
vdf = spdfcdfread(cdfdir,'Variables',{'vdf'});
pas2rtn_arr = spdfcdfread(cdfdir,'Variables',{'PAS_to_RTN'});

%% interpolate mag data
br_smooth = smoothdata(B_rtn(:,1),'movmean',501);
bt_smooth = smoothdata(B_rtn(:,2),'movmean',501);
bn_smooth = smoothdata(B_rtn(:,3),'movmean',501);
br_interp = interp1(epochmag,B_rtn(:,1),epoch);
bt_interp = interp1(epochmag,B_rtn(:,2),epoch);
bn_interp = interp1(epochmag,B_rtn(:,3),epoch);

babs_smooth = smoothdata(vecnorm(B_rtn,2,2),'movmean',501);
babs_interp = interp1(epochmag,vecnorm(B_rtn,2,2),epoch);
vA_interp = babs_interp*1e-9./sqrt(4e-7*pi*N_rtn*1e6*m_p)/1000;
%% plot vdf in RTN frame(ttime to select time)
sub_plot = find(epoch >= time_beg & epoch <= time_end);
epoch_plot = epoch(sub_plot);
walpha_arr = zeros(8,numel(sub_plot));
walpha_maxf_arr = zeros(8,numel(sub_plot));
tcount = 0;
is_lin_or_log_interp = 1;
%%
ttime = sub_plot(280);
    close all
    tcount = tcount+1;
    epoch_tmp = epoch(ttime);
    %% from pas to rtn coordinates
    pas2rtn = pas2rtn_arr(:,:,ttime);
    speedbin = sqrt(energybin*eV*2/m_p)/1000; %km/s
    temp_vdf = double((vdf(:,:,:,ttime)));
    [vsph,elesph,azisph] = meshgrid(speedbin,elebin,azibin);
    vx = -cosd(elesph).*vsph.*cosd(azisph);
    vy = -cosd(elesph).*vsph.*sind(azisph);
    vz = sind(elesph).*vsph;
    vr = vx*pas2rtn(1,1)+vy*pas2rtn(1,2)+vz*pas2rtn(1,3);
    vt = vx*pas2rtn(2,1)+vy*pas2rtn(2,2)+vz*pas2rtn(2,3);
    vn = vx*pas2rtn(3,1)+vy*pas2rtn(3,2)+vz*pas2rtn(3,3);
    % find max VDF in RTN coordinates
    maxVDF = max(temp_vdf,[],'all');
    submax = find(temp_vdf==maxVDF);
    vrmax_vdf = vr(submax); 
    vtmax_vdf = vt(submax); 
    vnmax_vdf = vn(submax);
    %get vdf>0
    is_vdf = temp_vdf>0;
     
    %calc interpolator
    vdf0 = temp_vdf(is_vdf); vr0 = vr(is_vdf); vt0 = vt(is_vdf); vn0 = vn(is_vdf);
    finterp = scatteredInterpolant(vr0,vt0,vn0,vdf0);
    finterp.ExtrapolationMethod =  'none';
    %interp to RTN mesh
    ngridr = 41;  ngridt = 17;  ngridn = 17;
    vrmin = 300;  vrmax = 700; dvr = (vrmax-vrmin)/(ngridr-1);
    vtmin = -200;  vtmax = 200; dvt = (vtmax-vtmin)/(ngridt-1);% better to include vy = 0
    vnmin = -200;  vnmax = 200; dvn = (vnmax-vnmin)/(ngridn-1);% better to include vz = 0
    vrn = linspace(vrmin,vrmax,ngridr);
    vtn = linspace(vtmin,vtmax,ngridt);
    vnn = linspace(vnmin,vnmax,ngridn);
    [vrnn,vtnn,vnnn]=meshgrid(vrn,vtn,vnn);
    VDFrtn = finterp(vrnn,vtnn,vnnn);
    %calc interpolator
    vdf0 = temp_vdf(is_vdf); vr0 = vr(is_vdf)-vrmax_vdf; vt0 = vt(is_vdf)-vtmax_vdf; vn0 = vn(is_vdf)-vnmax_vdf;
    finterp = scatteredInterpolant(vr0,vt0,vn0,vdf0);
    finterp.ExtrapolationMethod =  'none';
    %construct new mesh (vxn.vyn.vzn)for interp
    ngridx = 61;  ngridy = 61;  ngridz = 61;
    vxmin = -300;  vxmax = 300; dvx = (vxmax-vxmin)/(ngridx-1);
    vymin = -300;  vymax = 300; dvy = (vymax-vymin)/(ngridy-1);% better to include vy = 0
    vzmin = -300;  vzmax = 300; dvz = (vzmax-vzmin)/(ngridz-1);% better to include vz = 0
    vxn = linspace(vxmin,vxmax,ngridx);
    vyn = linspace(vymin,vymax,ngridy);
    vzn = linspace(vzmin,vzmax,ngridz);
    [vxnn,vynn,vznn]=meshgrid(vxn,vyn,vzn);
    %calc coordinates (vxnn,vynn,vznn) in rtn frame (vrnn,vtnn,vnnn)
    eb = [br_interp(ttime),bt_interp(ttime),bn_interp(ttime)];
    eb = eb/norm(eb);
    eperp2 = cross(eb,[1,0,0]);
    eperp2 = eperp2/norm(eperp2);
    eperp1 = cross(eperp2,eb);
    eperp1 = eperp1/norm(eperp1);
    vrnn = eb(1)*vxnn+eperp1(1)*vynn+eperp2(1)*vznn;
    vtnn = eb(2)*vxnn+eperp1(2)*vynn+eperp2(2)*vznn;
    vnnn = eb(3)*vxnn+eperp1(3)*vynn+eperp2(3)*vznn;
    VDFMFA = finterp(vrnn,vtnn,vnnn);    
%% plot VDF in RTN frame and MFA frame
figure;
b_plot_size = 100;
subplot(2,2,1)
contourf(vrn,vtn,log10(VDFrtn(:,:,(ngridn+1)/2)))
colormap jet
hold on
plot([vrmax_vdf vrmax_vdf+eb(1)*b_plot_size],[vtmax_vdf vtmax_vdf+eb(2)*b_plot_size],'k','LineWidth',2)
c=colorbar;
subplot(2,2,3)
contourf(vrn,vnn,log10(squeeze(VDFrtn((ngridt+1)/2,:,:))'))
colormap jet
subplot(2,2,2)
contourf(vxn,vyn,log10(VDFMFA(:,:,(ngridz+1)/2)))

%% overview plot
figure
subplot(5,1,1)
yyaxis left
plot(epochmag_plot,BLMN(1,:))
ylabel('B_L (nT)')
datetick('x','HH:MM')
yyaxis right
plot(epochmoms_plot,VLMN_moms(1,:))
ylabel('V_L (km/s)')
subplot(5,1,2)
yyaxis left
plot(epochmag_plot,BLMN(2,:))
ylabel('B_M (nT)')
datetick('x','HH:MM')
yyaxis right
plot(epochmoms_plot,VLMN_moms(2,:))
ylabel('V_M (km/s)')
subplot(5,1,3)
yyaxis left
plot(epochmag_plot,BLMN(3,:))
ylabel('B_N (nT)')
datetick('x','HH:MM')
yyaxis right
plot(epochmoms_plot,VLMN_moms(3,:))
ylabel('V_N (km/s)')
subplot(5,1,4)
plot(epochmag_plot,vecnorm(BLMN,2,1))
ylabel('|B| (nT)')
datetick('x','HH:MM')
subplot(5,1,5)
yyaxis left
plot(epochmoms_plot,N_rtn_plot)
ylabel('N (cm^{-3})')
datetick('x','HH:MM')
yyaxis right
plot(epochmoms_plot,T_rtn_plot)
ylabel('T (eV)')