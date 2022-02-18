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
for ttime = sub_plot(1):sub_plot(180)
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
    if is_lin_or_log_interp==1
        finterp = scatteredInterpolant(reshape(vr,[],1,1),reshape(vt,[],1,1),reshape(vn,[],1,1),reshape(temp_vdf,[],1,1));
    end
    if is_lin_or_log_interp==2
        finterp = scatteredInterpolant(reshape(vr,[],1,1),reshape(vt,[],1,1),reshape(vn,[],1,1),reshape(log10(temp_vdf),[],1,1));
    end
    finterp.ExtrapolationMethod =  'none';
    ngridr = 21;  ngridt = 21;  ngridn = 21;
    vrmin = 300;  vrmax = 700; dvr = (vrmax-vrmin)/(ngridr-1);
    vtmin = -200;  vtmax = 200; dvt = (vtmax-vtmin)/(ngridt-1);% better to include vy = 0
    vnmin = -200;  vnmax = 200; dvn = (vnmax-vnmin)/(ngridn-1);% better to include vz = 0
    vrn = linspace(vrmin,vrmax,ngridr);
    vtn = linspace(vtmin,vtmax,ngridt);
    vnn = linspace(vnmin,vnmax,ngridn);
    [vrnn,vtnn,vnnn]=meshgrid(vrn,vtn,vnn);
    VDFrtn = finterp(vrnn,vtnn,vnnn);
    if is_lin_or_log_interp==2
        VDFrtn = 10.^VDFrtn;
    end
    VDFrtn_write = VDFrtn;
    VDFrtn_write(isnan(VDFrtn_write)) = 0.0;
    VDFrtn(VDFrtn<1.e-12)=nan;
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
    temp_nflux = 2*VDFrtn/m_p/m_p.*(1000^2*0.5*m_p*(vrnn.^2+vtnn.^2+vnnn.^2)*eV)/100^4;
    vabsnn = sqrt(vrnn.^2+vtnn.^2+vnnn.^2);
    VDFrtn_v2 = smooth3(VDFrtn,'gaussian',3);
    sub_cal = VDFrtn>1e-15; % core + beam ?
     vr_cal = vrnn(sub_cal);
     vt_cal = vtnn(sub_cal);
     vn_cal = vnnn(sub_cal);
     vdf_cal = VDFrtn(sub_cal); 
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
    [lmax,mmax,nmax]=ind2sub([ngridt,ngridr,ngridn],find(VDFrtn==max(VDFrtn,[],'all')));
   % [VDFMFA,vxn,vxnn,vyn,vynn,vzn,vznn] = VDF_rtn2mfa(temp_vdf,vr,vt,vn,vrmax_vdf,vtmax_vdf,vnmax_vdf,eb,eperp1,eperp2);
%% interp to MFA frame using raw data
% finterp = scatteredInterpolant(reshape(vr,[],1,1),reshape(vt,[],1,1),reshape(vn,[],1,1),reshape(temp_vdf,[],1,1));
temp_int = temp_vdf > 0 & vr<540;
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
 [VDFpMFA,vxnp,vxnnp,vynp,vynnp,vznp,vznnp] = VDF_rtn2mfa(temp_vdf(vr<540),vr(vr<540),vt,vn,vrmax_vdf,vtmax_vdf,vnmax_vdf,eb,eperp1,eperp2);
VDFmfa(isnan(VDFmfa)|VDFmfa<1e-12)=0;
 v2d = trapz(vperp2n,VDFmfa,3);
v1d = trapz(vperp1n,v2d,1);
vpara_1d = v1d;
v1d = trapz(vparan,v2d,2);
vperp1_1d = v1d;
% plot(vparan,vpara_1d);
% hold on
% plot(vperp1n,vperp1_1d);
vdf_p_1d_para(:,tcount)=vpara_1d;
vdf_p_1d_perp(:,tcount)=vperp_1d;
end
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
    %% Subtract proton from VDF
    ngridr2 = 51;  ngridt2 = 21;  ngridn2 = 21;
    vrmin2 = 300;  vrmax2 = 700; dvr = (vrmax2-vrmin2)/(ngridr2-1);
    vtmin2 = -300;  vtmax2 = 300; dvt = (vtmax2-vtmin2)/(ngridt2-1);% better to include vy = 0
    vnmin2 = -300;  vnmax2 = 300; dvn = (vnmax2-vnmin2)/(ngridn2-1);% better to include vz = 0
    vrn2 = linspace(vrmin2,vrmax2,ngridr2);
    vtn2 = linspace(vtmin2,vtmax2,ngridt2);
    vnn2 = linspace(vnmin2,vnmax2,ngridn2);
    [vrnn2,vtnn2,vnnn2]=meshgrid(vrn2,vtn2,vnn2);
    VDFalpha = finterp(vrnn2,vtnn2,vnnn2);
    VDFalpha(VDFalpha<0|isnan(VDFalpha))=0;

    [maxf,maxn]=max(VDFalpha,[],'all','linear'); [indt,indr,indn]=ind2sub([ngridt2 ngridr2 ngridn2],maxn);
    %VDFalpha = VDFalpha(:,1:34,:);
    Ur_max = vrn2(indr); Ut_max = vtn2(indt); Un_max = vnn2(indn);
    figa = figure('visible','on');
    slice(vrn2,vtn2,vnn2,log10(VDFalpha),[550],[0],[0])
    xlabel('V_R(km/s)')
    ylabel('V_T(km/s)')
    zlabel('V_N(km/s)')
    colormap jet
    colorbar
    title(['f_a' datestr(epoch(ttime))])
%    savefig(figa,[savedir '\alpha_VDF_' datestr(epoch(ttime),'yyyymmdd_HHMMSS.fff') '.fig'])
    [wapara,waperp1,waperp2,na_mom,vra_mom,vta_mom,vna_mom,Q1] = calc_wpara_wperp(VDFalpha,vtn2,vrn2,vnn2,eb,eperp1,eperp2);
    walpha_arr(:,tcount)= [wapara,waperp1,waperp2,na_mom*1000,vra_mom,vta_mom,vna_mom,Q1];
     [wapara_max,waperp1_max,waperp2_max,na_mom,Q2] = calc_wpara_wperp_from_Uxyz(VDFalpha,vtn2,vrn2,vnn2,Ur_max,Ut_max,Un_max,eb,eperp1,eperp2);
walpha_maxf_arr(:,tcount)= [wapara_max,waperp1_max,waperp2_max,na_mom*1000,Ur_max,Ut_max,Un_max,Q2];

% beta0 = [0.01 0.0003 5 5 100];
% beta2 = nlinfit(vxn(v1d>0),log10(v1d(v1d>0)),@LogMaxwell_1D_fit,beta);
% figure;
% plot(vxn(v1d>0),log10(v1d(v1d>0)),'.')
% hold on
% plot(vxn(v1d>0),LogMaxwell_1D_fit(beta2,vxn(v1d>0)))
%% plot VDF slice || B0    
 [paramax,perp1max,perp2max]=ind2sub([61,41,41],find(VDFmfa==max(VDFmfa,[],'all')));
    is_lin_log = 1
    f5 = figure('visible','on');
    f5.Position =  [281 160.2000 1300 600];
    subplot(2,3,1)
    vrnplot = vrn;
    vtnplot = vtn;
    if is_lin_log == 1
        VDFplot = squeeze(log10(VDFrtn(:,:,nmax)));
        title_str = 'log_{10}VDF';
        crange = [-13 -8];
    end
    if is_lin_log == 2
        VDFplot = squeeze(VDFrtn(:,:,nmax));
        title_str = 'VDF';
        crange = 10.^[-13 -8];
    end
    contourf(vrnplot,vtnplot,VDFplot,[-8,-8.5,-9,-9.5,-10,-10.5,-11.5,-12.5,-13])
    grid on
   
    hold on
    p=plot([vrn(mmax),vrn(mmax)+eb(1)*300],...
         [vtn(lmax),vtn(lmax)+eb(2)*300],'k','LineWidth',2.5);
    xlabel('V_{R} (km/s)')
    ylabel('V_{T} (km/s)')
    %xlim([vrmin vrmax])
    shading flat
    colorbar
    colormap jet
    caxis(crange)
    title([title_str])

    subplot(2,3,4)
    vrnplot = vrn;
    vnnplot = vnn;
    if is_lin_log == 1
        VDFplot = squeeze(log10(VDFrtn(lmax,:,:))).';
    end
    if is_lin_log == 2
        VDFplot = squeeze(VDFrtn(lmax,:,:)).';
    end
    contourf(vrnplot,vnnplot,VDFplot,[-8,-8.5,-9,-9.5,-10,-10.5,-11.5,-12.5,-13])
    grid on
    hold on
    plot([vrn(mmax),vrn(mmax)+eb(1)*300],...
         [vnn(nmax),vnn(nmax)+eb(3)*300],'k','LineWidth',2.5)
    xlabel('V_{R} (km/s)')
    ylabel('V_{N} (km/s)')
   % xlim([vrmin vrmax])
    shading flat
    colorbar
    colormap jet
    caxis(crange)
    title([title_str])

    subplot(2,3,2)
    vrnplot = vparan;
    vtnplot = vperp1n;
    if is_lin_log == 1
        VDFplot = squeeze(log10(VDFmfa(:,:,perp2max)));
    end
    if is_lin_log == 2
        VDFplot = squeeze(VDFrtnfit(:,:,perp2max));
    end
    contourf(vrnplot,vtnplot,VDFplot,[-8,-8.5,-9,-9.5,-10,-10.5,-11.5,-12.5,-13])
    grid on
    xlabel('V_{//} (km/s)')
    ylabel('V_{\perp 1} (km/s)')
    shading flat
    colorbar
    colormap jet
    caxis(crange)
    title(['log_{10} VDF in MFA frame'])
    ylim([-200 200])
    xlim([-100 300])

    subplot(2,3,5)
    vrnplot = vparan;
    vnnplot = vperp2n;
    if is_lin_log == 1
        VDFplot = squeeze(log10(VDFmfa(perp1max,:,:))).';
    end
    if is_lin_log == 2
        VDFplot = squeeze(VDFrtnfit(lmax,:,:)).';
    end
    contourf(vrnplot,vnnplot,VDFplot,[-8,-8.5,-9,-9.5,-10,-10.5,-11.5,-12.5,-13])
    grid on
    xlabel('V_{//} (km/s)')
    ylabel('V_{\perp 2} (km/s)')
    shading flat
    colorbar
    colormap jet
    caxis(crange)
    title(['log_{10} VDF in MFA frame'])
    
    subplot(2,3,3)
    p1=plot(epochmag_plot,BLMN(1,:),'k');
    hold on;
    p2=plot(epochmag_plot,BLMN(2,:),'b');
    p3=plot(epochmag_plot,BLMN(3,:),'r');
    ylabel('B (nT)')
    xlabel(date_str)
    yl=ylim;
    xlim([mag_time_beg mag_time_end])
    p4=plot([epoch_tmp,epoch_tmp],yl,'--r','LineWidth',2);
    datetick('x','HH:MM:SS','keeplimits')
    legend([p1 p2 p3],'B_L','B_M','B_N','Locations','northeast');
    
    subplot(2,3,6)
    yyaxis left
    p1=plot(epochmoms_plot,N_rtn_plot,'Linewidth',1.5);
    ylabel('n (cm^{-3})')
    yyaxis right
    %p2=plot(epochmoms_plot,T_rtn_plot,'Linewidth',1.5);
    %ylabel('T (eV)')
    p2=plot(epochmoms_plot,v_A_arr/1000,'Linewidth',1.5);
    ylabel('v_A(km/s)')
    xlabel(date_str)
    yl=ylim;
    xlim([mag_time_beg mag_time_end])
    hold on;
    p3=plot([epoch_tmp,epoch_tmp],yl,'--r','LineWidth',2);
    datetick('x','HH:MM:SS','keeplimits')
    
    %saveas(f2,[savedir '\2D_VDF_rtn_fitting_' datestr(epoch(ttime),'yyyymmdd_HHMMSS.fff') '.png'])

%% plot VDF slice || B0    
    is_lin_log = 1
    f2 = figure('visible','on');
    f2.Position =  [281 160.2000 1300 600];
    subplot(2,3,1)
    vrnplot = vrn;
    vtnplot = vtn;
    if is_lin_log == 1
        VDFplot = squeeze(log10(VDFrtn(:,:,nmax)));
        title_str = 'log_{10}VDF';
        crange = [-13 -8];
    end
    if is_lin_log == 2
        VDFplot = squeeze(VDFrtn(:,:,nmax));
        title_str = 'VDF';
        crange = 10.^[-13 -8];
    end
    contourf(vrnplot,vtnplot,VDFplot,[-8,-8.5,-9,-9.5,-10,-10.5,-11.5,-12.5,-14.5])
    grid on
   
    hold on
    p=plot([vrn(mmax),vrn(mmax)+eb(1)*300],...
         [vtn(lmax),vtn(lmax)+eb(2)*300],'k','LineWidth',2.5);
    xlabel('V_{R} (km/s)')
    ylabel('V_{T} (km/s)')
    xlim([vrmin vrmax])
    shading flat
    colorbar
    colormap jet
    caxis(crange)
    title([title_str])

    subplot(2,3,4)
    vrnplot = vrn;
    vnnplot = vnn;
    if is_lin_log == 1
        VDFplot = squeeze(log10(VDFrtn(lmax,:,:))).';
    end
    if is_lin_log == 2
        VDFplot = squeeze(VDFrtn(lmax,:,:)).';
    end
    contourf(vrnplot,vnnplot,VDFplot,[-8,-8.5,-9,-9.5,-10,-10.5,-11.5,-12.5,-14.5])
    grid on
    hold on
    plot([vrn(mmax),vrn(mmax)+eb(1)*300],...
         [vnn(nmax),vnn(nmax)+eb(3)*300],'k','LineWidth',2.5)
    xlabel('V_{R} (km/s)')
    ylabel('V_{N} (km/s)')
    xlim([vrmin vrmax])
    shading flat
    colorbar
    colormap jet
    caxis(crange)
    title([title_str])

    subplot(2,3,2)
    vrnplot = vrn;
    vtnplot = vtn;
    if is_lin_log == 1
        VDFplot = squeeze(log10(VDFrtnfit(:,:,nmax)));
    end
    if is_lin_log == 2
        VDFplot = squeeze(VDFrtnfit(:,:,nmax));
    end
    contourf(vrnplot,vtnplot,VDFplot,[-8,-8.5,-9,-9.5,-10,-10.5,-11.5,-12.5,-14.5])
    grid on
    xlabel('V_{R} (km/s)')
    ylabel('V_{T} (km/s)')
    shading flat
    colorbar
    colormap jet
    caxis(crange)
    title(['【fit result】log_{10} VDF'])

    subplot(2,3,5)
    vrnplot = vrn;
    vnnplot = vnn;
    if is_lin_log == 1
        VDFplot = squeeze(log10(VDFrtnfit(lmax,:,:))).';
    end
    if is_lin_log == 2
        VDFplot = squeeze(VDFrtnfit(lmax,:,:)).';
    end
    contourf(vrnplot,vnnplot,VDFplot,[-8,-8.5,-9,-9.5,-10,-10.5,-11.5,-12.5,-14.5])
    grid on
    xlabel('V_{R} (km/s)')
    ylabel('V_{N} (km/s)')
    shading flat
    colorbar
    colormap jet
    caxis(crange)
    title(['【fit result】log_{10} VDF'])
    
    subplot(2,3,3)
    p1=plot(epochmag_plot,BLMN(1,:),'k');
    hold on;
    p2=plot(epochmag_plot,BLMN(2,:),'b');
    p3=plot(epochmag_plot,BLMN(3,:),'r');
    ylabel('B (nT)')
    xlabel(date_str)
    yl=ylim;
    xlim([mag_time_beg mag_time_end])
    p4=plot([epoch_tmp,epoch_tmp],yl,'--r','LineWidth',2);
    datetick('x','HH:MM:SS','keeplimits')
    legend([p1 p2 p3],'B_L','B_M','B_N','Locations','northeast');
    
    subplot(2,3,6)
    yyaxis left
    p1=plot(epochmoms_plot,N_rtn_plot,'Linewidth',1.5);
    ylabel('n (cm^{-3})')
    yyaxis right
    %p2=plot(epochmoms_plot,T_rtn_plot,'Linewidth',1.5);
    %ylabel('T (eV)')
    p2=plot(epochmoms_plot,v_A_arr/1000,'Linewidth',1.5);
    ylabel('v_A(km/s)')
    xlabel(date_str)
    yl=ylim;
    xlim([mag_time_beg mag_time_end])
    hold on;
    p3=plot([epoch_tmp,epoch_tmp],yl,'--r','LineWidth',2);
    datetick('x','HH:MM:SS','keeplimits')
    
    saveas(f2,[savedir '\2D_VDF_rtn_fitting_onlycb_' datestr(epoch(ttime),'yyyymmdd_HHMMSS.fff') '.png'])
    %% plot 1D VDF along B0
    f3 = figure('visible','on');
    f3.Position =  [281 160.2000 1200 401.8000];
    subplot(1,3,1)
    vrnplot = vrn;
    if is_lin_log == 1
        dataplot = squeeze(log10(VDFrtn(lmax,:,nmax)));
        fitplot = squeeze(log10(VDFrtnfit(lmax,:,nmax)));
        fitcoreplot = squeeze(log10(VDFcorefit(lmax,:,nmax)));
        fitbeamplot = squeeze(log10(VDFbeamfit(lmax,:,nmax)));
        if is_cb_or_cba == 2
            fitalphaplot = squeeze(log10(VDFalphafit(lmax,:,nmax)));
        end
    end
    if is_lin_log == 2
        dataplot = squeeze(VDFrtn(lmax,:,nmax));
        fitplot = squeeze(VDFrtnfit(lmax,:,nmax));
        fitcoreplot = squeeze(VDFcorefit(lmax,:,nmax));
        fitbeamplot = squeeze(VDFbeamfit(lmax,:,nmax));
        if is_cb_or_cba == 2
            fitalphaplot = squeeze(VDFalphafit(lmax,:,nmax));
        end
    end
    if is_cb_or_cba == 1
        plot(vrnplot,dataplot,'+m','MarkerSize',15)
        hold on
        plot(vrnplot,fitcoreplot,'r',...
             vrnplot,fitbeamplot,'b',...
             vrnplot,fitplot,'--k','LineWidth',1.5)
    end
    if is_cb_or_cba == 2
        plot(vrnplot(vrnplot>530)/sqrt(2),dataplot(vrnplot>530),'+m','MarkerSize',15)
        hold on
       % plot(vrnplot,fitcoreplot,'r',...
        %     vrnplot,fitbeamplot,'b',...
         %    vrnplot,fitalphaplot,'g',...
          %   vrnplot,fitplot,'--k','LineWidth',1.5)
    end
    if is_lin_log==1
        ylim([-12 -6.5])
    end
    xlim([350 500])
    xlabel('V_{//} (km/s)')
    if is_cb_or_cba==1
        legend('data','core','beam','total','Location','northeast')
    end
    if is_cb_or_cba==2
        legend('data','core','beam','alpha','total','Location','northeast')
    end
    ylabel('log_{10}VDF')

    subplot(1,3,2)
    vtnplot = vtn;
    if is_lin_log == 1
        dataplot = squeeze(log10(VDFrtn(:,mmax,nmax)));
        fitplot = squeeze(log10(VDFrtnfit(:,mmax,nmax)));
        fitcoreplot = squeeze(log10(VDFcorefit(:,mmax,nmax)));
        fitbeamplot = squeeze(log10(VDFbeamfit(:,mmax,nmax)));
        if is_cb_or_cba == 2
            fitalphaplot = squeeze(log10(VDFalphafit(:,mmax,nmax)));
        end
    end
    if is_lin_log == 2
        dataplot = squeeze(VDFrtn(:,mmax,nmax));
        fitplot = squeeze(VDFxyzfit(:,mmax,nmax));
        fitcoreplot = squeeze(VDFcorefit(:,mmax,nmax));
        fitbeamplot = squeeze(VDFbeamfit(:,mmax,nmax));
        if is_cb_or_cba == 2
            fitalphaplot = squeeze(VDFalphafit(:,mmax,nmax));
        end
    end
    if is_cb_or_cba == 1
        plot(vtnplot,dataplot,'+m','MarkerSize',15)
        hold on
        plot(vtnplot,fitcoreplot,'r',...
             vtnplot,fitbeamplot,'b',...
             vtnplot,fitplot,'--k','LineWidth',1.5)
    end
    if is_cb_or_cba == 2
        plot(vtnplot,dataplot,'+m','MarkerSize',15)
        hold on
        plot(vtnplot,fitcoreplot,'r',...
             vtnplot,fitbeamplot,'b',...
             vtnplot,fitalphaplot,'g',...
             vtnplot,fitplot,'--k','LineWidth',1.5)
    end
    if is_lin_log ==1
        ylim([-12 -6.5])
    end
    xlabel('V_{T} (km/s)')
%     if is_cb_or_cba==1
%         legend('data','core','beam','total','Location','northwest')
%     end
%     if is_cb_or_cba==2
%         legend('data','core','beam','alpha','total','Location','northwest')
%     end
    ylabel('log_{10}VDF')


    subplot(1,3,3)
    vnnplot = vnn;
    if is_lin_log == 1
        dataplot = squeeze(log10(VDFrtn(lmax,mmax,:)));
        fitplot = squeeze(log10(VDFrtnfit(lmax,mmax,:)));
        fitcoreplot = squeeze(log10(VDFcorefit(lmax,mmax,:)));
        fitbeamplot = squeeze(log10(VDFbeamfit(lmax,mmax,:)));
    end
    if is_lin_log == 2
        dataplot = squeeze(VDFrtn(lmax,mmax,:));
        fitplot = squeeze(VDFrtnfit(lmax,mmax,:));
        fitcoreplot = squeeze(VDFcorefit(lmax,mmax,:));
        fitbeamplot = squeeze(VDFbeamfit(lmax,mmax,:));
    end
    if is_cb_or_cba == 1
        plot(vnnplot,dataplot,'+m','MarkerSize',15)
        hold on
        plot(vnnplot,fitcoreplot,'r',...
             vnnplot,fitbeamplot,'b',...
             vnnplot,fitplot,'--k','LineWidth',1.5)
    end
    if is_cb_or_cba == 2
        plot(vnnplot,dataplot,'+m','MarkerSize',15)
        hold on
        plot(vnnplot,fitcoreplot,'r',...
             vnnplot,fitbeamplot,'b',...
             vnnplot,fitalphaplot,'g',...
             vnnplot,fitplot,'--k','LineWidth',1.5)
    end
    if is_lin_log ==1
        ylim([-12 -6.5])
    end
    xlabel('V_{N} (km/s)')
%     if is_cb_or_cba==1
%         legend('data','core','beam','total','Location','northwest')
%     end
%     if is_cb_or_cba==2
%         legend('data','core','beam','alpha','total','Location','northwest')
%     end
    ylabel('log_{10}VDF')
     saveas(f3,[savedir       '\1D_VDF_fitting_onlycb_' datestr(epoch(ttime),'yyyymmdd_HHMMSS.fff') '.png'])
     %%
%      f4 = figure('visible','off');
%         x1 = vr_cal; x2 = vt_cal; x3 = vn_cal;
%         VDFcorefitori = 1.e-3*c(1)/(c(4)^2*c(5)*pi^1.5).*exp(-((x1-c(10)).*eb(1)+(x2-c(11)).*eb(2)+(x3-c(12)).*eb(3)).^2/c(5)^2).*exp(-((x1-c(10)).*eperp1(1)+(x2-c(11)).*eperp1(2)+(x3-c(12)).*eperp1(3)).^2/c(4)^2).*exp(-((x1-c(10)).*eperp2(1)+(x2-c(11)).*eperp2(2)+(x3-c(12)).*eperp2(3)).^2/c(4)^2);
%         VDFbeamfitori = 1.e-3*c(2)/(c(6)^2*c(7)*pi^1.5).*exp(-((x1-c(13)).*eb(1)+(x2-c(14)).*eb(2)+(x3-c(15)).*eb(3)).^2/c(7)^2).*exp(-((x1-c(13)).*eperp1(1)+(x2-c(14)).*eperp1(2)+(x3-c(15)).*eperp1(3)).^2/c(6)^2).*exp(-((x1-c(13)).*eperp2(1)+(x2-c(14)).*eperp2(2)+(x3-c(15)).*eperp2(3)).^2/c(6)^2);
%         VDFalphafitori = 1.e-3*c(3)/(c(8)^2*c(9)*pi^1.5).*exp(-((x1-c(16)).*eb(1)+(x2-c(17)).*eb(2)+(x3-c(18)).*eb(3)).^2/c(9)^2).*exp(-((x1-c(16)).*eperp1(1)+(x2-c(17)).*eperp1(2)+(x3-c(18)).*eperp1(3)).^2/c(8)^2).*exp(-((x1-c(16)).*eperp2(1)+(x2-c(17)).*eperp2(2)+(x3-c(18)).*eperp2(3)).^2/c(8)^2);
%         VDFrtnfitori = VDFcorefitori+VDFbeamfitori+VDFalphafitori;
%      subplot(1,3,1)
% semilogy(vr_cal,vdf_cal,'.')
% hold on
% semilogy(x1,VDFrtnfitori,'.')
% xlabel('V_R')
% ylim([10e-12,10e-8])
% legend('original','fit')
% subplot(1,3,2)
% semilogy(vt_cal,vdf_cal,'.')
% hold on
% semilogy(x2,VDFrtnfitori,'.')
% xlabel('V_T')
% ylim([10e-12,10e-8])
% subplot(1,3,3)
% semilogy(vn_cal,vdf_cal,'.')
% hold on
% semilogy(x3,VDFrtnfitori,'.')
% xlabel('V_N')
% ylim([10e-12,10e-8])
% f4.Position = [125         343        1024         420];
%  saveas(f4,[savedir       '\1D_VDF_original_' datestr(epoch(ttime),'yyyymmdd_HHMMSS.fff') '.png'])
end
%%
load(['coeff_alpha.mat'])
f1 = figure;
plot(epoch(sub_plot),walpha_arr(1:3,:))
datetick('x','HH:MM','keeplimits')

yyaxis right
plot(epoch(sub_plot),br_interp(sub_plot),'k','LineWidth',2)

legend('w_{a//}','w_{a\perp 1}','w_{a\perp 2}','B_R')

f2 = figure;
plot(epoch(sub_plot),(walpha_arr(5:7,:)/sqrt(2)-coeff_arr(:,10:12)')./vA_interp(sub_plot)')
datetick('x','HH:MM','keeplimits')
title('U_\alpha/1.414-U_c')
yyaxis right
plot(epoch(sub_plot),br_interp(sub_plot),'k','LineWidth',2)

legend('U_{dR}','U_{dT}','U_{dN}','B_R')

f3 = figure;
vd_mom = walpha_arr(5:7,:)/sqrt(2)-coeff_arr(:,10:12)';
Bvec_interp = [br_interp bt_interp bn_interp]';
thetavdB = acosd(dot(vd_mom,Bvec_interp(:,sub_plot),1)./vecnorm(vd_mom,2,1)./vecnorm(Bvec_interp(:,sub_plot),2,1));
plot(epoch(sub_plot),thetavdB,'LineWidth',2)
datetick('x','HH:MM','keeplimits')
title('\theta(U_d,B)')
yyaxis right
plot(epoch(sub_plot),br_interp(sub_plot),'k','LineWidth',2)

f4 = figure;
plot(epoch(sub_plot),coeff_arr(:,4:5),'LineWidth',1)
yyaxis right
plot(epoch(sub_plot),br_interp(sub_plot),'k','LineWidth',2)
datetick('x','HH:MM','keeplimits')
title('w_c')
legend('w_{c\perp}','w_{c//}','B_R')
%%
load([savedir  '\fit_coeff.mat'])
subplot(4,1,1)
plot(epoch(sub_plot),br_interp(sub_plot))
hold on
plot(epoch(sub_plot),bt_interp(sub_plot))
plot(epoch(sub_plot),bn_interp(sub_plot))
legend('B_R','B_T','B_N')
datetick('x','HH:MM','keeplimits')
subplot(4,1,2)
plot(epoch_arr,coeff_arr(:,1:3))
hold on
plot(epoch(sub_plot),N_rtn(sub_plot))
datetick('x','HH:MM','keeplimits')
legend('n_c','n_b','n_a','n_p')
subplot(4,1,3)
plot(epoch_arr,coeff_arr(:,4:5))
hold on
plot(epoch_arr,smoothdata(coeff_arr(:,8:9)))
datetick('x','HH:MM','keeplimits')
legend('T_cpara','T_cperp','T_apara','T_aperp')
subplot(4,1,4)
dV = sqrt((coeff_arr(:,16)/2-coeff_arr(:,10)).^2+(coeff_arr(:,17)/2-coeff_arr(:,11)).^2+(coeff_arr(:,1)/2-coeff_arr(:,12).^2));
plot(epoch_arr,dV)
datetick('x','HH:MM','keeplimits')


Vcore = coeff_arr(:,10:12);
%% save fitting parameters
%     filesave = ['/Users/psr/work/离子回旋波统计/Data/' event_type '/(high_res)ICW_fitting_parameters_' ...
%         datestr(time_beg,'yyyymmddHHMMSS') '-' datestr(time_end,'yyyymmddHHMMSS') '.dat'];
%     save(filesave,'coeff_arr','-ascii');
filesave = ['E:\work\离子回旋波统计\Data\' event_type '\ICW_fitting_parameters_' ...
        datestr(time_beg,'yyyymmddHHMMSS') '-' datestr(time_end,'yyyymmddHHMMSS') '.dat'];
    save(filesave,'coeff_arr','-ascii');