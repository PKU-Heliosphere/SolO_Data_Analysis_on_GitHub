% PAS_VDF_RTN_ga_fitting
%% for meeting
figure
p1=plot(epochmag_plot,BLMN(1,:),'k');
    hold on;
    p2=plot(epochmag_plot,BLMN(2,:),'b');
    p3=plot(epochmag_plot,BLMN(3,:),'r');
    ylabel('B (nT)')
    xlabel(date_str)
    yl=ylim;
    xlim([mag_time_beg mag_time_end])
    plot([epoch(sub_plot(80)),epoch(sub_plot(80))],yl,'--r','LineWidth',2);
    plot([epoch(sub_plot(220)),epoch(sub_plot(220))],yl,'--r','LineWidth',2);
     plot([epoch(sub_plot(250)),epoch(sub_plot(250))],yl,'--r','LineWidth',2);
      plot([epoch(sub_plot(380)),epoch(sub_plot(380))],yl,'--r','LineWidth',2);
    datetick('x','HH:MM:SS','keeplimits')
    legend([p1 p2 p3],'B_L','B_M','B_N','Locations','northeast');
%

%% This script is used to plot SolO PAS high resolution and fitting in RTN coordinates
%close all;
%clc; clear all; %#ok<CLALL>
event_type = 'magnetic reconnection'
onesec = datenum(2010,1,1,1,1,1)-datenum(2010,1,1,1,1,0);
m_p = 1.673e-27;%kg
eV = 1.602e-19;%J
date_str = '20201014';
hour_beg = 22; min_beg = 45; sec_beg = 00;
hour_end = 23; min_end = 10; sec_end = 00;
year_str = date_str(1:4); mon_str = date_str(5:6); day_str = date_str(7:8);
year = str2num(year_str); month = str2num(mon_str); day = str2num(day_str);
time_beg = datenum(year,month,day,hour_beg,min_beg,sec_beg);
time_end = datenum(year,month,day,hour_end,min_end,sec_end);
mag_hour_beg = 22; mag_min_beg = 45; mag_sec_beg = 00;
mag_hour_end = 23; mag_min_end = 10; mag_sec_end = 00;
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
int_LMN = [1:100 460:510];
br_plot = BRp; bt_plot = BTp; bn_plot = BNp;
epochmag_plot = subepochp;
[BL,BM,BN] = calc_mag_MVA(br_plot(int_LMN),bt_plot(int_LMN),bn_plot(int_LMN));
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
wproton_arr = zeros(8,numel(sub_plot));
walpha_arr = zeros(8,numel(sub_plot));
walpha_maxf_arr = zeros(8,numel(sub_plot));
epara_arr = zeros(3,numel(sub_plot));
eperp1_arr = zeros(3,numel(sub_plot));
eperp2_arr = zeros(3,numel(sub_plot));
vdf_p_1d_para = zeros(61,numel(sub_plot));
vdf_p_1d_perp = zeros(41,numel(sub_plot));
vdf_a_1d_para = zeros(61,numel(sub_plot));
vdf_a_1d_perp = zeros(41,numel(sub_plot));
tcount = 0;
is_lin_or_log_interp = 1;
%%
p_or_alpha = 'p';
for ttime = sub_plot(1):sub_plot(510)
    close all
    tcount = ttime-sub_plot(1)+1;
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
%     subnan = find(temp_vdf==0);
%     temp_vdf(temp_vdf==0)=-1;
    %% calc moments
%     int_dV = -vsph.^2.*cosd(elesph).*1e6;
%     int_phi = trapz(deg2rad(azibin),temp_vdf.*int_dV,3);
%     int_theta = trapz(deg2rad(elebin),int_phi,1);
%     int_v = trapz(speedbin*1000,int_theta);
    
    %% interpolate to a regular grid in rtn coordinates
%     if is_lin_or_log_interp==1
%         finterp = scatteredInterpolant(reshape(vr,[],1,1),reshape(vt,[],1,1),reshape(vn,[],1,1),reshape(temp_vdf,[],1,1));
%     end
%     if is_lin_or_log_interp==2
%         finterp = scatteredInterpolant(reshape(vr,[],1,1),reshape(vt,[],1,1),reshape(vn,[],1,1),reshape(log10(temp_vdf),[],1,1));
%     end
%     finterp.ExtrapolationMethod =  'none';
%     ngridr = 21;  ngridt = 21;  ngridn = 21;
%     vrmin = 300;  vrmax = 700; dvr = (vrmax-vrmin)/(ngridr-1);
%     vtmin = -200;  vtmax = 200; dvt = (vtmax-vtmin)/(ngridt-1);% better to include vy = 0
%     vnmin = -200;  vnmax = 200; dvn = (vnmax-vnmin)/(ngridn-1);% better to include vz = 0
%     vrn = linspace(vrmin,vrmax,ngridr);
%     vtn = linspace(vtmin,vtmax,ngridt);
%     vnn = linspace(vnmin,vnmax,ngridn);
%     [vrnn,vtnn,vnnn]=meshgrid(vrn,vtn,vnn);
%     VDFrtn = finterp(vrnn,vtnn,vnnn);
%     if is_lin_or_log_interp==2
%         VDFrtn = 10.^VDFrtn;
%     end
%     VDFrtn_write = VDFrtn;
%     VDFrtn_write(isnan(VDFrtn_write)) = 0.0;
%     VDFrtn(VDFrtn<1.e-12)=nan;
%     VDFrtn = smooth3(VDFrtn,'gaussian',2);
  
%     ngridr = 101;  ngridt = 51;  ngridn = 51;
%     vrmin = 200;  vrmax = 600; dvr = (vrmax-vrmin)/(ngridr-1);
%     vtmin = -200;  vtmax = 200; dvt = (vtmax-vtmin)/(ngridt-1);% better to include vy = 0
%     vnmin = -200;  vnmax = 200; dvn = (vnmax-vnmin)/(ngridn-1);% better to include vz = 0
%     vrn = linspace(vrmin,vrmax,ngridr);
%     vtn = linspace(vtmin,vtmax,ngridt);
%     vnn = linspace(vnmin,vnmax,ngridn);
%     [vrnn,vtnn,vnnn]=meshgrid(vrn,vtn,vnn);
%     vr_data = reshape(vr,[],1,1);
%     vt_data = reshape(vt,[],1,1);
%     vn_data = reshape(vn,[],1,1);
%     VDF_data = reshape(log10(temp_vdf),[],1,1);
%     vr_data = vr_data(~isinf(VDF_data));
%     vt_data = vt_data(~isinf(VDF_data));
%     vn_data = vn_data(~isinf(VDF_data));
%     VDF_data = VDF_data(~isinf(VDF_data));
%     VDFrtn = griddata(vr_data,vt_data,vn_data,VDF_data,...
%                   vrnn,vtnn,vnnn,'natural');
%     VDFrtn = 10.^VDFrtn;
%     VDFrtn_write = VDFrtn;
%     VDFrtn_write(isnan(VDFrtn_write)) = 0.0;
%     VDFrtn(VDFrtn==0)=nan;
    
%     slice(vrnn,vtnn,vnnn,log10(VDFrtn),[0 116],0,0)
    %% get mean magnetic field vector
    eb = [br_interp(ttime),bt_interp(ttime),bn_interp(ttime)];
    eb = eb/norm(eb); epara_arr(:,tcount)=eb;
    eperp2 = cross(eb,[1,0,0]);
    eperp2 = eperp2/norm(eperp2); eperp2_arr(:,tcount)=eperp2;
    eperp1 = cross(eperp2,eb);
    eperp1 = eperp1/norm(eperp1); eperp1_arr(:,tcount)=eperp1;
%     temp_nflux = 2*VDFrtn/m_p/m_p.*(1000^2*0.5*m_p*(vrnn.^2+vtnn.^2+vnnn.^2)*eV)/100^4;
%     vabsnn = sqrt(vrnn.^2+vtnn.^2+vnnn.^2);
%     VDFrtn_v2 = smooth3(VDFrtn,'gaussian',3);
%     sub_cal = VDFrtn>1e-15; % core + beam ?
%      vr_cal = vrnn(sub_cal);
%      vt_cal = vtnn(sub_cal);
%      vn_cal = vnnn(sub_cal);
%      vdf_cal = VDFrtn(sub_cal); 
     %vdf_cal(vr_cal>600)=0;
%      
% vr_cal0 = reshape(vr,[],1,1);
% vt_cal0 = reshape(vt,[],1,1);
% vn_cal0 = reshape(vn,[],1,1);
% vdf_cal0 = reshape(temp_vdf,[],1,1);
% temp_nflux = reshape(2*temp_vdf/m_p/m_p.*(1000^2*0.5*m_p*(vr.^2+vn.^2+vt.^2)*eV)/100^4,[],1,1);
% %fitflag =  temp_nflux>10 & vr_cal0 < 600 & vr_cal0 > 330; 
% fitflag = vdf_cal0>10^(-10);
% vr_cal = vr_cal0(fitflag); vt_cal = vt_cal0(fitflag); vn_cal = vn_cal0(fitflag); vdf_cal=vdf_cal0(fitflag);
%    [lmax,mmax,nmax]=ind2sub([ngridt,ngridr,ngridn],find(VDFrtn==max(VDFrtn,[],'all')));
   % [VDFMFA,vxn,vxnn,vyn,vynn,vzn,vznn] = VDF_rtn2mfa(temp_vdf,vr,vt,vn,vrmax_vdf,vtmax_vdf,vnmax_vdf,eb,eperp1,eperp2);
%% interp to MFA frame using raw data
% finterp = scatteredInterpolant(reshape(vr,[],1,1),reshape(vt,[],1,1),reshape(vn,[],1,1),reshape(temp_vdf,[],1,1));
vabs = sqrt(vr.^2+vt.^2+vn.^2);
temp_1d_vdf = (eflux(ttime,35:50));
peakloc=find(imregionalmin(temp_1d_vdf))+35-1;
posi_index = find(temp_1d_vdf>0)+35-1;
if numel(peakloc)==1
vpa = speedbin(peakloc(end));
vaa = speedbin(posi_index(1));
end
if numel(peakloc)>1
    vpa=speedbin(peakloc(end)); vaa = speedbin(peakloc(1));
end
if p_or_alpha == 'a'
temp_int = temp_vdf > 0 & vabs>=vpa & vabs<=vaa;
else
    temp_int = temp_vdf > 0 & vabs<vpa;
end
vr_temp =vr(temp_int); vt_temp = vt(temp_int); vn_temp = vn(temp_int);
RTN_2_MFA = [eb;eperp1;eperp2];
vmfa_temp = RTN_2_MFA*[vr_temp vt_temp vn_temp]';
[maxvdf,maxindex]=max(temp_vdf(temp_int));
finterp_mfa = scatteredInterpolant(vmfa_temp(1,:)'-vmfa_temp(1,maxindex),vmfa_temp(2,:)'-vmfa_temp(2,maxindex),vmfa_temp(3,:)'-vmfa_temp(3,maxindex),temp_vdf(temp_int));

% interp to B para perp frame    
    npara = 61;  nperp1 = 41;  nperp2 = 41;
    vparamin = -300;  vparamax = 300; dvpara = (vparamax-vparamin)/(npara-1);
    vperp1min = -200;  vperp1max = 200; dvperp1 = (vperp1max-vperp1min)/(nperp1-1);% better to include vy = 0
    vperp2min = -200;  vperp2max = 200; dvperp2 = (vperp2max-vperp2min)/(nperp2-1);% better to include vz = 0
    vparan = linspace(vparamin,vparamax,npara);
    vperp1n = linspace(vperp1min,vperp1max,nperp1);
    vperp2n = linspace(vperp2min,vperp2max,nperp2);
    [vparann,vperp1nn,vperp2nn]=meshgrid(vparan,vperp1n,vperp2n);
    finterp_mfa.ExtrapolationMethod='none';
    VDFmfa = finterp_mfa(vparann,vperp1nn,vperp2nn);
    VDFmfa(VDFmfa<1e-12)=0;

%  plot mfa frame
%     xplot = vparan;
%     yplot = vperp1n;
%     fplot = log10(VDFmfa(:,:,21));
%     contourf(xplot,yplot,fplot,[-8,-8.5,-9,-9.5,-10,-10.5,-11.5,-12.5,-14.5])

    % calc 1D VDF 
% [VDFpMFA,vxnp,vxnnp,vynp,vynnp,vznp,vznnp] = VDF_rtn2mfa(temp_vdf(vr<540),vr(vr<540),vt,vn,vrmax_vdf,vtmax_vdf,vnmax_vdf,eb,eperp1,eperp2);
if p_or_alpha == 'a'
    vparan = vparan/sqrt(2);
    vperp1n = vperp1n/sqrt(2);
    vperp2n = vperp2n/sqrt(2);
end
VDFmfa(isnan(VDFmfa)|VDFmfa<1e-12)=0;
v2d = trapz(vperp2n,VDFmfa,3);
v1d = trapz(vperp1n,v2d,1);
vpara_1d = v1d;
v1d = trapz(vparan,v2d,2);
vperp1_1d = v1d;
% plot(vparan,vpara_1d);
% hold on
% plot(vperp1n,vperp1_1d);
if p_or_alpha == 'a'
 [wapara,waperp1,waperp2,na_mom,vra_mom,vta_mom,vna_mom,Q1] = calc_wpara_wperp(VDFmfa,vperp1n,vparan,vperp2n,eb,eperp1,eperp2,4,'T');
else
    [wapara,waperp1,waperp2,na_mom,vra_mom,vta_mom,vna_mom,Q1] = calc_wpara_wperp(VDFmfa,vperp1n,vparan,vperp2n,eb,eperp1,eperp2,1,'T');
end
v_max = vmfa_temp(:,maxindex);
v_max = v_max+[vra_mom,vta_mom,vna_mom]';
v_max_rtn = inv(RTN_2_MFA)*v_max;
if p_or_alpha == 'a' 
v_max_rtn = v_max_rtn/sqrt(2);
vdf_a_1d_para(:,tcount)=vpara_1d;
vdf_a_1d_perp(:,tcount)=vperp1_1d;
walpha_arr(:,tcount)= [wapara,waperp1,waperp2,na_mom*1000,v_max_rtn(1),v_max_rtn(2),v_max_rtn(3),Q1];
else
vdf_p_1d_para(:,tcount)=vpara_1d;
vdf_p_1d_perp(:,tcount)=vperp1_1d;
wproton_arr(:,tcount)= [wapara,waperp1,waperp2,na_mom*1000,v_max_rtn(1),v_max_rtn(2),v_max_rtn(3),Q1];
end
end

%% Plot 1d vdf
t = tiledlayout(7,1,'TileSpacing','none');
nexttile
plotrange = 1:510;
plotx = subepochp(plotrange);
pcolor(plotx,vparan,log10(vdf_p_1d_para(:,plotrange)))
shading flat
ylim([-200 200])
caxis([-7 -3])
colorbar
datetick('x','HH:MM')
xticklabels([])
ylabel('V_{p//}')
nexttile
pcolor(plotx,vperp1n,log10(vdf_p_1d_perp(:,plotrange)))
shading flat
datetick('x','HH:MM')
xticklabels([])
ylim([-200 200])
caxis([-7 -3])
%colorbar
ylabel('V_{p\perp}')

nexttile
pcolor(plotx,vparan/sqrt(2),log10(vdf_a_1d_para(:,plotrange)))
shading flat
datetick('x','HH:MM')
xticklabels([])
ylim([-150 150])
caxis([-7 -4.5])
ylabel('V_{\alpha //}')

nexttile
pcolor(plotx,vperp1n/sqrt(2),log10(vdf_a_1d_perp(:,plotrange)))
shading flat
ylim([-150 150])
datetick('x','HH:MM')
xticklabels([])
colormap jet
c = colorbar;
%c.Layout.Tile = 'east';
caxis([-7 -4.5])
ylabel('V_{\alpha\perp}')





nexttile
yyaxis left
plot(subepochp,wproton_arr(4,:))
ylabel('n_p (cm^{-3})')
yyaxis right
plot(subepochp,walpha_arr(4,:))
ylabel('n_\alpha (cm^{-3})')
ylim([0 0.4])
datetick('x','HH:MM')
xticklabels([])
nexttile
vp_LMN = [BL;BM;BN]*wproton_arr(5:7,:);
va_LMN = [BL;BM;BN]*walpha_arr(5:7,:);
plot(subepochp,vp_LMN(1,:))
hold on
plot(subepochp,va_LMN(1,:))
ylabel('V_L (km/s)')
legend('V_{pL}','V_{\alpha L}')
datetick('x','HH:MM')
xticklabels([])
nexttile
yyaxis left
l1 = plot(subepochp,wproton_arr(1,:),'LineWidth',1);
hold on
l2 = plot(subepochp,(wproton_arr(2,:)+wproton_arr(3,:))/2,'LineWidth',1);
ylim([0 25])
ylabel('T_p (eV)')
yyaxis right
l3 = plot(subepochp,walpha_arr(1,:),'LineWidth',1);
l4 = plot(subepochp,(walpha_arr(2,:)+walpha_arr(3,:))/2,'LineWidth',1);
ylim([0 70])
ylabel('T_\alpha (eV)')

datetick('x','HH:MM')
hold on
line1_x = [subepochp(104) subepochp(104)];
line2_x = [subepochp(180) subepochp(180)];
line3_x = [subepochp(426) subepochp(426)];
line4_x = [subepochp(291) subepochp(291)];
line5_x = [subepochp(246) subepochp(246)];
line6_x = [subepochp(267) subepochp(267)];
line_y = [0 70*7];
line(line1_x,line_y,'LineStyle','--','LineWidth',2,'Color','#D95319','Clipping','off')
line(line2_x,line_y,'LineStyle','--','LineWidth',2,'Color',	'#7E2F8E','Clipping','off')
line(line3_x,line_y,'LineStyle','--','LineWidth',2,'Color','#D95319','Clipping','off')
line(line4_x,line_y,'LineStyle','--','LineWidth',2,'Color',	'#7E2F8E','Clipping','off')
line(line5_x,line_y,'LineStyle','--','LineWidth',2,'Color',		'#77AC30','Clipping','off')
line(line6_x,line_y,'LineStyle','--','LineWidth',2,'Color',		'#77AC30','Clipping','off')
legend([l1 l2 l3 l4],{'T_{p//}','T_{p\perp}','T_{\alpha //}','T_{\alpha \perp}'})
%%

%% GA method fitting
%     x11=vr_cal;
%     x22=vt_cal;
%     x33=vn_cal;
%     yy=vdf_cal;
%     tic
%     save tempdata.mat x11 x22 x33 yy eb eperp1 eperp2;
%     is_cb_or_cba = 1;
%     % c(12) =
%     % [nc,nb,vcthp,vcthz,vbthp,vbthz,vcswpara,vcswperp1,vcswperp2,vbswpara,vbswperp1,vbswperp2]
%     if is_cb_or_cba == 1
%         LB = [2,0.01,...
%               10,10,10,10,...
%               300,-100,-100,...
%               300,-100,-100];
%         UB = [20,10,...
%               100,100,100,100,...
%               500,100,100,... 
%               500,150,100];
%         ObjectiveFunction = @ga_cb_RTN_fitness;
%         nvars = 12; % number of varibles
%         ConstraintFunction = []; % constraints
%         rng default; % for reproducibality ?
%         options=optimoptions('ga','FunctionTolerance',1.e-5);
%         [coeff,fval]=ga(ObjectiveFunction,nvars,...
%             [-1 2 0 0 0 0 0 0 0 0 0 0],0,[],[],LB,UB,ConstraintFunction,options);
%         c=coeff;
%         VDFcorefit = 1.e-3*c(1)/(c(3)^2*c(4)*pi^1.5)*exp(-(vrnn-c(7)).^2/c(4)^2).*exp(-(vtnn-c(8)).^2/c(3)^2).*exp(-(vnnn-c(9)).^2/c(3)^2);
%         VDFbeamfit = 1.e-3*c(2)/(c(5)^2*c(6)*pi^1.5)*exp(-(vrnn-c(10)).^2/c(6)^2).*exp(-(vtnn-c(11)).^2/c(5)^2).*exp(-(vnnn-c(12)).^2/c(5)^2);
%         VDFrtnfit = 1.e-3*c(1)/(c(3)^2*c(4)*pi^1.5)*exp(-(vrnn-c(7)).^2/c(4)^2).*exp(-(vtnn-c(8)).^2/c(3)^2).*exp(-(vnnn-c(9)).^2/c(3)^2)+...
%             1.e-3*c(2)/(c(5)^2*c(6)*pi^1.5)*exp(-(vrnn-c(10)).^2/c(6)^2).*exp(-(vtnn-c(11)).^2/c(5)^2).*exp(-(vnnn-c(12)).^2/c(5)^2);
% %        VDFrtnfit(VDFxyzfit<1.e-15)=nan;
%     end
%     % c(18) =
%     % [nc,nb,na,...
%     %  vcthp,vcthz,vbthp,vbthz,vathp,vathz,...
%     %  vcswpara,vcswperp1,vcswperp2,...
%     %  vbswpara,vbswperp1,vbswperp2,...
%     %  vcaswpara,vaswperp1,vaswperp2]
%     if is_cb_or_cba == 2
%         LB = [3,0.1,0.1,...
%               10,10,10,10,10,10,...
%               350,-100,-100,...
%              300,-100,-100,...
%               500,-100,-100];
%         UB = [20,10,5,...
%               80,80,100,100,100,100,...
%               500,100,100,...
%               550,100,100,...
%               630,100,100];
%         ObjectiveFunction = @ga_cba_RTN_fitness;
%         nvars = 18; % number of varibles
%         ConstraintFunction = []; % constraints
%         rng default; % for reproducibality ?
%        % options=optimoptions('ga','FunctionTolerance',1.e-8);
%         [coeff,fval]=ga(ObjectiveFunction,nvars,...
%             [],[],[],[],LB,UB,ConstraintFunction);
%         
%         c=coeff;
%         x1 = vrnn; x2 = vtnn; x3 = vnnn;
%         VDFcorefit = 1.e-3*c(1)/(c(4)^2*c(5)*pi^1.5).*exp(-((x1-c(10)).*eb(1)+(x2-c(11)).*eb(2)+(x3-c(12)).*eb(3)).^2/c(5)^2).*exp(-((x1-c(10)).*eperp1(1)+(x2-c(11)).*eperp1(2)+(x3-c(12)).*eperp1(3)).^2/c(4)^2).*exp(-((x1-c(10)).*eperp2(1)+(x2-c(11)).*eperp2(2)+(x3-c(12)).*eperp2(3)).^2/c(4)^2);
%         VDFbeamfit = 1.e-3*c(2)/(c(6)^2*c(7)*pi^1.5).*exp(-((x1-c(13)).*eb(1)+(x2-c(14)).*eb(2)+(x3-c(15)).*eb(3)).^2/c(7)^2).*exp(-((x1-c(13)).*eperp1(1)+(x2-c(14)).*eperp1(2)+(x3-c(15)).*eperp1(3)).^2/c(6)^2).*exp(-((x1-c(13)).*eperp2(1)+(x2-c(14)).*eperp2(2)+(x3-c(15)).*eperp2(3)).^2/c(6)^2);
%         VDFalphafit = 1.e-3*c(3)/(c(8)^2*c(9)*pi^1.5).*exp(-((x1-c(16)).*eb(1)+(x2-c(17)).*eb(2)+(x3-c(18)).*eb(3)).^2/c(9)^2).*exp(-((x1-c(16)).*eperp1(1)+(x2-c(17)).*eperp1(2)+(x3-c(18)).*eperp1(3)).^2/c(8)^2).*exp(-((x1-c(16)).*eperp2(1)+(x2-c(17)).*eperp2(2)+(x3-c(18)).*eperp2(3)).^2/c(8)^2);
%         VDFrtnfit = VDFcorefit+VDFbeamfit+VDFalphafit;
%         VDFrtnfit(VDFrtnfit<1.e-19)=nan;
%     end
%         if is_cb_or_cba == 3
%         LB = [3,0.2,0.2,0.2,...
%               10,10,10,10,10,10,10,10,...
%               350,-100,-100,...
%              300,-100,-100,...
%               500,-100,-100,...
%               500,-100,-100];
%         UB = [20,10,2,2,...
%               80,80,100,100,100,100,100,100,...
%               500,100,100,...
%               550,100,100,...
%               630,100,100,...
%               700,100,100];
%         ObjectiveFunction = @ga_cbab_RTN_fitness;
%         nvars = 24; % number of varibles
%         ConstraintFunction = []; % constraints
%         rng default; % for reproducibality ?
%         [coeff,fval]=ga(ObjectiveFunction,nvars,...
%             [],[],[],[],LB,UB,ConstraintFunction);
%         options=optimoptions(@ga,'FunctionTolerance',1.e-30,'PopulationSize',100000);
%         c=coeff;
%         x1 = vrnn; x2 = vtnn; x3 = vnnn;
%         VDFcorefit = 1.e-3*c(1)/(c(5)^2*c(6)*pi^1.5).*exp(-((x1-c(13)).*eb(1)+(x2-c(14)).*eb(2)+(x3-c(15)).*eb(3)).^2/c(5)^2).*exp(-((x1-c(13)).*eperp1(1)+(x2-c(14)).*eperp1(2)+(x3-c(15)).*eperp1(3)).^2/c(6)^2).*exp(-((x1-c(13)).*eperp2(1)+(x2-c(14)).*eperp2(2)+(x3-c(15)).*eperp2(3)).^2/c(6)^2);
%         VDFbeamfit = 1.e-3*c(2)/(c(7)^2*c(8)*pi^1.5).*exp(-((x1-c(16)).*eb(1)+(x2-c(17)).*eb(2)+(x3-c(18)).*eb(3)).^2/c(7)^2).*exp(-((x1-c(16)).*eperp1(1)+(x2-c(17)).*eperp1(2)+(x3-c(18)).*eperp1(3)).^2/c(8)^2).*exp(-((x1-c(16)).*eperp2(1)+(x2-c(17)).*eperp2(2)+(x3-c(18)).*eperp2(3)).^2/c(8)^2);
%         VDFalphafit =1.e-3*c(3)/(c(9)^2*c(10)*pi^1.5).*exp(-((x1-c(19)).*eb(1)+(x2-c(20)).*eb(2)+(x3-c(21)).*eb(3)).^2/c(9)^2).*exp(-((x1-c(19)).*eperp1(1)+(x2-c(20)).*eperp1(2)+(x3-c(21)).*eperp1(3)).^2/c(10)^2).*exp(-((x1-c(19)).*eperp2(1)+(x2-c(20)).*eperp2(2)+(x3-c(21)).*eperp2(3)).^2/c(10)^2);
%         VDFalphabeamfit = 1.e-3*c(4)/(c(11)^2*c(12)*pi^1.5).*exp(-((x1-c(22)).*eb(1)+(x2-c(23)).*eb(2)+(x3-c(24)).*eb(3)).^2/c(11)^2).*exp(-((x1-c(22)).*eperp1(1)+(x2-c(23)).*eperp1(2)+(x3-c(24)).*eperp1(3)).^2/c(12)^2).*exp(-((x1-c(22)).*eperp2(1)+(x2-c(23)).*eperp2(2)+(x3-c(24)).*eperp2(3)).^2/c(12)^2);
%         VDFrtnfit = VDFcorefit+VDFbeamfit+VDFalphafit+VDFalphabeamfit;
%         VDFrtnfit(VDFrtnfit<1.e-19)=nan;
%         end
%     toc
%     if is_cb_or_cba == 1
%         coeff_arr(tcount,1:12) = coeff;
%         coeff_arr(tcount,13:15) = [vrn(mmax),vtn(lmax),vnn(nmax)];
%         coeff_arr(tcount,16:18) = [br_interp(ttime),bt_interp(ttime),bn_interp(ttime)];
%     end
%     if is_cb_or_cba == 2
%         coeff_arr(tcount,1:18) = coeff;
%         coeff_arr(tcount,19:21) = [vrn(mmax),vtn(lmax),vnn(nmax)];
%         coeff_arr(tcount,22:24) = [br_interp(ttime),bt_interp(ttime),bn_interp(ttime)];
%     end
