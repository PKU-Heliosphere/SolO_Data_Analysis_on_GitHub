%% This script is used to calculate the HT frame using SolO data
mu_0 = 4e-7*pi;
m_p = 1.673e-27;
k_B = 1.380649e-23;

date_str = '20201014';
year = str2num(date_str(1:4)); month = str2num(date_str(5:6)); day = str2num(date_str(7:8));
hour_beg = 22; min_beg = 45; sec_beg = 0; time_beg = datenum(year,month,day,hour_beg,min_beg,sec_beg);
hour_end = 23; min_end = 10; sec_end = 0; time_end = datenum(year,month,day,hour_end,min_end,sec_end);

solodir = 'D:\SolOData\';
easname = 'solo_l2_swa-eas2-nm3d-psd_20201014t000035-20201014t235855_v01.cdf';

cdfdir = [solodir 'solo_l2_swa-pas-vdf_' date_str '_v01.cdf'];
magdir = [solodir 'solo_l2_mag-rtn-normal_' date_str '_v02.cdf'];
if ~isfile(magdir)
    magdir = [solodir 'solo_l2_mag-rtn-normal_' date_str '_v01.cdf'];
end
pmomdir = [solodir 'solo_l2_swa-pas-grnd-mom_' date_str '_v01.cdf'];
pasefluxdir = [solodir 'solo_l2_swa-pas-eflux_' date_str '_v01.cdf'];
easpadefluxdir = [solodir easname];
onesec = datenum(2010,1,1,1,1,1)-datenum(2010,1,1,1,1,0);
m_p = 1.673e-27;%kg
eV = 1.602e-19;%J


epochmag = spdfcdfread(magdir,'Variables',{'EPOCH'});
B_rtn = double(spdfcdfread(magdir,'Variables',{'B_RTN'}));
epochp = spdfcdfread(pmomdir,'Variables',{'Epoch'});
np = spdfcdfread(pmomdir,'Variables',{'N'});
vp = spdfcdfread(pmomdir,'Variables',{'V_RTN'});
Pp = spdfcdfread(pmomdir,'Variables',{'P_RTN'});
Tp = spdfcdfread(pmomdir,'Variables',{'T'});

%% load 1D ion eflux
efluxepoch = spdfcdfread(pasefluxdir,'Variables',{'Epoch'});
eflux = spdfcdfread(pasefluxdir,'Variables',{'eflux'});

%% subinterval and interpolate
subint = epochp >= time_beg & epochp <= time_end; subepochp = epochp(subint);
BRp = interp1(epochmag,B_rtn(:,1),subepochp);
BTp = interp1(epochmag,B_rtn(:,2),subepochp);
BNp = interp1(epochmag,B_rtn(:,3),subepochp);
% subint = epochmag >= time_beg & epochmag <= time_end; subepochp = epochmag(subint);
% BRp = B_rtn(subint,1); BTp = B_rtn(subint,2); BNp = B_rtn(subint,3);

% vp_mag = interp1(epochp,vp,epochmag);
% vp = vp_mag;
 Babsp = sqrt(BRp.^2+BTp.^2+BNp.^2);
vA_arr = Babsp*1e-9./sqrt(4e-7*pi*np(subint)*1e6*m_p);
cs_arr = sqrt(5*Tp(subint)*1.6e-19/3/m_p);
%%
plot(epochmag,B_rtn)
hold on
plot(epochmag,sqrt(sum(B_rtn.^2,2)))
datetick('x','HH:MM:SS')
legend('B_R','B_T','B_N','|B|')
title('B (nT) 0.8AU')

%% Plot V & B
load('Data\Reconnection_Paper_Plotting\alpha_moment.mat')
load('Data\Reconnection_Paper_Plotting\proton_moment.mat')
%n0 = mean(np(subint));
B0 = mean(Babsp);
f1=figure(1)
t = tiledlayout(8,1,'TileSpacing','none');
nexttile % |B|
plot(subepochp,sqrt(BRp.^2+BTp.^2+BNp.^2),'k-','LineWidth',1);
ylim([5.5 10])
ylabel('|B| (nT)')
datetick('x','HH:MM')
xticklabels([])
title(date_str)
nexttile %BRTN
thetaBR = acosd(BRp./Babsp);
l1 = plot(subepochp,BRp);
subint_index = find(subint);
subint_beg = subint_index(1); subint_end = subint_index(end);
energypeak = zeros(numel(subint_index),2);
for i=subint_beg:subint_end
   [pks,locs]=findpeaks(log10(eflux(i,:)),'MinPeakDistance',9,'SortStr','descend');
   if i == 661+subint_beg-1
       [pks,locs]=findpeaks(log10(eflux(i,:)),'MinPeakDistance',8,'SortStr','descend');
   end
   energypeak(i-subint_beg+1,1)=energybin(locs(1));
   energypeak(i-subint_beg+1,2)=energybin(locs(2));
end

plot(subepochp,[BRp BTp BNp])
ylabel('B_{RTN} (nT)')
%plot([epochp(17097) epochp(17097)],[300 400],'r--','LineWidth',1)
% plot([epochp(17090) epochp(17090)],[300 400],'r--','LineWidth',1)
% plot([epochp(17110) epochp(17110)],[300 400],'r--','LineWidth',1)
% plot([epochp(17260) epochp(17260)],[300 400],'r--','LineWidth',1)
% plot([epochp(17280) epochp(17280)],[300 400],'r--','LineWidth',1)
% plot([epochp(17300) epochp(17300)],[300 400],'r--','LineWidth',1)

datetick('x','HH:MM')
legend({'B_R','B_T','B_N'},'Location','eastoutside')
xticklabels([])

nexttile

plot(subepochp,vp(subint,1)-mean(vp(subint,1),1))
hold on
plot(subepochp,vp(subint,2:3))
ylabel('V_i (km/s)')
legend({'V_R - 411','V_T','V_N'},'Location','eastoutside')
datetick('x','HH:MM')
xticklabels([])


nexttile

plot(subepochp,BLMN)
ylabel('B_{LMN} (nT)')
datetick('x','HH:MM')
legend({'B_L','B_M','B_N'},'Location','eastoutside')
ylim([-11 11])
xticklabels([])
% subplot(614)
% 
% pcolor(tt(tt>=time_beg & tt <= time_end),f(fwave),ds_sigma_m(:,tt>=time_beg & tt <= time_end))
% title('\sigma_{mTN} (red = 1, blue = -1)')
% shading flat
% set(gca,'yscale','log')
% datetick('x','HH:MM')
% ylabel('f (Hz)')
% yticks([0.01 0.1 0.5])
% yticklabels({'0.01','0.1','1'})
% colormap jet
% yyaxis right
% plot(subepochp,thetaBR,'k')
% ylabel('\theta_{BR}')
% ylim([0 180])
% yticks([0 90 180])
% yticklabels({'0','90','180'})
nexttile

plot(subepochp,VLMN-mean(VLMN,2))
ylabel('V_{LMN} (nT)')
datetick('x','HH:MM')
xticklabels([])
ylim([-50 44])
legend({'V_L + 270','V_T - 154','V_N + 270'},'Location','eastoutside')
nexttile
yyaxis left

plot(subepochp,wproton_arr(4,:))
%plot(subepochp,np(subint))
ylabel('n_i (cm^{-3})')
yyaxis right
%plot(subepochp,Tp(subint))
plot(subepochp,wproton_arr(1,:),'r')
hold on

plot(subepochp,0.5*(wproton_arr(2,:)+wproton_arr(3,:)),'k-')
ylabel('T_i (eV)')
datetick('x','HH:MM')
xticklabels([])
legend('n_i','T_{i//}','T_{i\perp}','Location','eastoutside')
nexttile

speedplot = speedbin <800 & speedbin > 200;
pcolor(subepochp,(energybin(speedplot)),log10(eflux(subint,speedplot)'))
hold on
[maxnum,imaxmax]=max(eflux(subint,:)');
plot(subepochp,(energypeak(:,1)),'w','LineWidth',1)
plot(subepochp,(energypeak(:,2)),'k','LineWidth',1)
shading flat
%title('1D ion Eflux (cm^{-2}s^{-1}eV/eV)')
ylabel('Energy (eV/e)')

datetick('x','HH:MM')
xticklabels([])
colormap jet
c1 = colorbar;
c1.Label.String = 'E flux (cm^{-2}s^{-1}eV/eV)';
c1.Ticks = [6 7 8 9];
c1.TickLabels = {'10^6','10^7','10^8','10^9'};
%c1.LineWidth = 1; 
c1.TickLength=0.05;
nexttile
easint = easepoch>=time_beg & easepoch<=time_end;
pcolor(easepoch(easint),linspace(0,180,n_PAbin),log10(squeeze(sum(strahlPA(easint,1:32,:),2)))');
xticks([0 90 180])
xticklabels({'0','90','180'})
hold on
%ylim([0 180])
%[maxnum,imaxmax]=max(eflux(subint,:)');
%plot(subepochp,speedbin(imaxmax),'w')
shading flat
c2 = colorbar;
c2.Label.String = 'PSD (s^3km^{-6})';
c2.Ticks = [3 4];
c2.TickLabels = {'10^3','10^4'};
c2.TickLength = 0.05;
%colormap(h2,'jet')
%title('93 eV PA PSD')
ylabel('Pitch Angle')
datetick('x','HH:MM')

line1_x = [subepochp(180) subepochp(180)];
line2_x = [subepochp(291) subepochp(291)];
line3_x = [subepochp(426) subepochp(426)];
line4_x = [subepochp(206) subepochp(206)];
line_y = [0 180*8];
line(line1_x,line_y,'LineStyle','--','LineWidth',1,'Color','#D95319','Clipping','off')
line(line2_x,line_y,'LineStyle','--','LineWidth',1,'Color',	'#7E2F8E','Clipping','off')
line(line3_x,line_y,'LineStyle','--','LineWidth',1,'Color','#D95319','Clipping','off')
line(line4_x,line_y,'LineStyle','--','LineWidth',1,'Color',	'#7E2F8E','Clipping','off')
%% Figure 2 LMN
int_LMN = [1:100 460:510];
br_plot = BRp; bt_plot = BTp; bn_plot = BNp;
epochmag_plot = subepochp;
[BL,BM,BN] = calc_mag_MVA(br_plot(int_LMN),bt_plot(int_LMN),bn_plot(int_LMN));
BLMN = [BL;BM;BN]*[br_plot bt_plot bn_plot]';
figure;
subplot(3,1,1)
plot(epochmag_plot,[br_plot bt_plot bn_plot])
hold on
plot(epochmag_plot,vecnorm([br_plot bt_plot bn_plot],2,2))
datetick('x','HH:MM')
legend('B_R','B_T','B_N','|B|')
subplot(3,1,2)
plot(epochmag_plot,BLMN)
hold on
plot(epochmag_plot,vecnorm([br_plot bt_plot bn_plot],2,2))
datetick('x','HH:MM')
legend('B_L','B_M','B_N','|B|')
subplot(3,1,3)
VLMN = [BL;BM;BN]*vp(subint,:)';
plot(epochmag_plot,VLMN-mean(VLMN,2))
datetick('x','HH:MM')
%%
plot([br_plot bt_plot bn_plot])
hold on
plot(vecnorm([br_plot bt_plot bn_plot],2,2))
%%
%HTint = 225:261; %RD1
HTint = 261:280; %RD2
%HTint = 85:122;
limint = [-70 70];
Bint = BLMN(:,HTint); Vint = VLMN(:,HTint);
BHTintx = BLMN(1,HTint); BHTinty=BLMN(2,HTint); BHTintz=BLMN(3,HTint);
VHTintx = VLMN(1,HTint); VHTinty=VLMN(2,HTint); VHTintz=VLMN(3,HTint);
vHT = calc_HT_frame(BHTintx,BHTinty,BHTintz,VHTintx,VHTinty,VHTintz);
VinHTframe = Vint-vHT;
Eadv = cross(-Vint,Bint,1);
EHT = cross(repmat(vHT,1,size(Bint,2)),Bint,1);
fit_quality = fit(Eadv(1,:)',EHT(1,:)','poly1');
vA_vec = BLMN*1e-9./sqrt(4e-7*pi*np(subint)'*1e6*m_p)/1000;
vA_HT = vA_vec(:,HTint);
subplot(2,2,1)
plot(vA_HT(1,:),Vint(1,:)-vHT(1),'o');
fit_WT1 = fit(vA_HT(1,:)',Vint(1,:)'-vHT(1),'poly1');
fit_WT2 = fit(vA_HT(2,:)',Vint(2,:)'-vHT(2),'poly1');
fit_WT3 = fit(vA_HT(3,:)',Vint(3,:)'-vHT(3),'poly1');
hold on
plot(vA_HT(2,:),Vint(2,:)-vHT(2),'o');
plot(vA_HT(3,:),Vint(3,:)-vHT(3),'o');
legend({['CC_{VL,VAL}=' num2str(fit_WT1.p1,'%.2f')],['CC_{VM,VAM}=' num2str(fit_WT2.p1,'%.2f')],['CC_{VN,VAN}=' num2str(fit_WT3.p1,'%.2f')]})
xticklabels
xlabel('V_A (km/s)','FontSize',18)
ylabel('V-V_{HT} (km/s)','FontSize',18)
%plot([0 50],-[0 50],'k')
%plot(vA_HT/1000,Vint(3,:)-vHT(3),'o');
title('Walen Test')
xlim(limint); ylim(limint)
axis square
subplot(2,2,2)
limint = [-2500 2500];
plot(Eadv',EHT','o')
title('HT frame Quality')
axis square
xlabel('E=-V\times B','FontSize',18)
ylabel('E_{HT}=-V_{HT}\times B','FontSize',18)
fit_Q1 = fit(Eadv(1,:)',EHT(1,:)','poly1');
fit_Q2 = fit(Eadv(2,:)',EHT(2,:)','poly1');
fit_Q3 = fit(Eadv(3,:)',EHT(3,:)','poly1');
legend({['CC_{EL,EHTL}=' num2str(fit_Q1.p1,'%.2f')],['CC_{EM,EHTM}=' num2str(fit_Q2.p1,'%.2f')],['CC_{EN,EHTN}=' num2str(fit_Q3.p1,'%.2f')]})
xlim(limint); ylim(limint)
axis square
subplot(2,2,[3 4])
yyaxis left
plot(Vint'-vHT')
ylabel('V_{HTframe} (km/s)')
yyaxis right
l2=plot(Bint');
xlabel('sample number')
ylabel('B (nT)')
legend(l2,'B_R','B_T','B_N')
%% determine RD normal direction
[BLRD1,BMRD1,BNRD1] = calc_mag_MVA(br_plot(HTint),bt_plot(HTint),bn_plot(HTint));
BLMNRD = [BLRD1;BMRD1;BNRD1]*Bint;
plot(BLMNRD')
%%
shock_index  = findchangepts(sqrt(BRp.^2+BTp.^2+BNp.^2),'MaxNumChanges',5,'Statistic','mean');
% shock_index = [104 182 206 266 291]
VRp = vp(subint,1); VTp = vp(subint,2); VNp = vp(subint,3); Np_subint = np(subint); Tp_subint = Tp(subint);
PB = (BRp.^2+BTp.^2+BNp.^2)*1e-18/2/mu_0;
Pkinetic = Np_subint*1e6.*Tp_subint*11605*k_B;

%int_before = 91:96; int_after = 108:113; %DC1
%int_before = 167:172; int_after = 189:194; %DC2
%int_before = 195:200; int_after = 210:4:230; %DC3
%int_before = 250:260; int_after = 273:283; %DC4
    int_before = 216:230; int_after = 231:243; %DC4
Bbefore = [BRp(int_before) BTp(int_before) BNp(int_before)]*1e-9;
Vbefore = [VRp(int_before) VTp(int_before) VNp(int_before)]*1e3;
nbefore = Np_subint(int_before)*1e6; Pbefore = Pkinetic(int_before);
Bafter = [BRp(int_after) BTp(int_after) BNp(int_after)]*1e-9;
Vafter = [VRp(int_after) VTp(int_after) VNp(int_after)]*1e3;
nafter = Np_subint(int_after)*1e6; Pafter = Pkinetic(int_after);
%[n_theta,n_phi,V_shock] = Shock_Fitting_2D_t_meshgrid(Bbefore,Vbefore,nbefore,Pbefore,Bafter,Vafter,nafter,Pafter);
%%
plotint = int_before(1):int_after(end);
B1x = Bbefore(:,1); B1y = Bbefore(:,2); B1z = Bbefore(:,3);
V1x = Vbefore(:,1); V1y = Vbefore(:,2); V1z = Vbefore(:,3);
B2x = Bafter(:,1); B2y = Bafter(:,2); B2z = Bafter(:,3);
V2x = Vafter(:,1); V2y = Vafter(:,2); V2z = Vafter(:,3);
n1 = nbefore; n2 = nafter; P1 = Pbefore; P2 = Pafter;
nx = cosd(n_theta); ny = sind(n_theta).*cosd(n_phi); nz = sind(n_theta).*sind(n_phi);
e_sfx = [nx,ny,nz]; e_sfy = cross([1,0,0],e_sfx); e_sfy=e_sfy/norm(e_sfy);
e_sfz = cross(e_sfx,e_sfy); e_sfz = e_sfz/norm(e_sfz);

B_sf = [e_sfx;e_sfy;e_sfz]*[BRp(plotint) BTp(plotint) BNp(plotint)]';
V_sf = [e_sfx;e_sfy;e_sfz]*[VRp(plotint) VTp(plotint) VNp(plotint)]';
Vs = V_shock;
massflux = Np_subint(plotint).*(V_sf(1,:)-Vs/1000)';

Bn1 = B1x.*nx+B1y.*ny+B1z.*nz; Bn1x = Bn1.*nx; Bn1y = Bn1.*ny; Bn1z = Bn1.*nz;
Bn2 = B2x.*nx+B2y.*ny+B2z.*nz; Bn2x = Bn2.*nx; Bn2y = Bn2.*ny; Bn2z = Bn2.*nz;
Bt1x = B1x-Bn1x; Bt1y = B1y-Bn1y; Bt1z = B1z-Bn1z; Bt1sqr = Bt1x.^2+Bt1y.^2+Bt1z.^2;
Bt2x = B2x-Bn2x; Bt2y = B2y-Bn2y; Bt2z = B2z-Bn2z; Bt2sqr = Bt2x.^2+Bt2y.^2+Bt2z.^2;
Vn1 = V1x.*nx+V1y.*ny+V1z.*nz; Vn1x = Vn1.*nx; Vn1y = Vn1.*ny; Vn1z = Vn1.*nz;
Vn2 = V2x.*nx+V2y.*ny+V2z.*nz; Vn2x = Vn2.*nx; Vn2y = Vn2.*ny; Vn2z = Vn2.*nz;
Vt1x = V1x-Vn1x; Vt1y = V1y-Vn1y; Vt1z = V1z-Vn1z; %Vt1 = sqrt(Vt1x.^2+Vt1y.^2+Vt1z.^2);
Vt2x = V2x-Vn2x; Vt2y = V2y-Vn2y; Vt2z = V2z-Vn2z; %Vt2 = sqrt(Vt2x.^2+Vt2y.^2+Vt2z.^2);
dGn = m_p.*(n1.*(Vn1-Vs)-n2.*(Vn2-Vs));
dSn= P1+Bt1sqr/2/mu_0+m_p.*n1.*(Vn1-Vs).^2-...
    (P2+Bt2sqr/2/mu_0+m_p.*n2.*(Vn2-Vs).^2);

% Plot in Shock frame
fshock = tiledlayout(6,1);
%plot |B|,n

nexttile

plotint = int_before(1):int_after(end);
plotx = subepochp(plotint);
ploty = BRp(plotint).^2+BTp(plotint).^2+BNp(plotint).^2;
plotBabs = sqrt(ploty);
ploty = plotBabs;
yyaxis left
plot(plotx,ploty)
hold on
plot(subepochp(int_before(end))*ones(1,2),[min(ploty) max(ploty)])
plot(subepochp(int_after(1))*ones(1,2),[min(ploty) max(ploty)])
ylabel('|B| (nT)','FontSize',16)
yyaxis right

ploty = Np_subint(plotint);
plot(plotx,ploty);
ylabel('n_i (cm^{-3})','FontSize',14)
datetick('x','HH:MM')
xticklabels([])
%plot Bn, Vn-Vs
nexttile

yyaxis left
ploty = [B_sf(1,:)];
plot(plotx,ploty)
ylabel('B_n (nT)','FontSize',16)
datetick('x','HH:MM')
yyaxis right
plot(plotx,(V_sf(1,:)-Vs/1000))
ylabel('V_n-V_s (km/s)','FontSize',12)
xticklabels([])
nexttile
%plot Bt,Vt
yyaxis left
ploty =sqrt(B_sf(2,:).^2+B_sf(3,:).^2);
plot(plotx,ploty)
ylabel('B_t (nT)','FontSize',16)
yyaxis right
ploty = sqrt(V_sf(2,:).^2+V_sf(3,:).^2);
plot(plotx,ploty)
ylabel('v_t (km/s)','FontSize',14)
datetick('x','HH:MM')
xticklabels([])
nexttile
%plot T,ptot
yyaxis left
ploty = Tp_subint(plotint);

plot(plotx,ploty)
ylabel('T (eV)','FontSize',16)
yyaxis right
ploty = PB(plotint)+Pkinetic(plotint);
plot(plotx,ploty)
ylabel('p_{tot}=P_i+P_B (Pa)','FontSize',12)
datetick('x','HH:MM')
xticklabels([])
nexttile
%plot thetaBn
thetaBn = acosd((BRp(plotint).*nx+BTp(plotint).*ny+BNp(plotint).*nz)./plotBabs);
ploty = thetaBn;
plot(plotx,ploty)
ylabel('\theta_{Bn}','FontSize',16)
datetick('x','HH:MM')
xticklabels([])
nexttile
%plot VA, Cf, Cs

vAplot =vA_arr(plotint)/1000; csplot = cs_arr(plotint)/1000;
cfplot = 0.5*((vAplot.^2+csplot.^2)+sqrt((vAplot.^2+csplot.^2).^2-4*csplot.^2.*vAplot.^2.*cosd(thetaBn).^2));
cslowplot = 0.5*((vAplot.^2+csplot.^2)-sqrt((vAplot.^2+csplot.^2).^2-4*csplot.^2.*vAplot.^2.*cosd(thetaBn).^2));
plot(plotx,abs(abs(V_sf(1,:)-Vs/1000)'./vAplot./cosd(thetaBn)))
hold on
plot(plotx,abs(V_sf(1,:)-Vs/1000)'./sqrt(cfplot))
plot(plotx,abs(V_sf(1,:)-Vs/1000)'./sqrt(cslowplot))
plot(plotx,ones(size(plotx)),'k')
%plot(plotx,abs(V_sf(1,:)-Vs/1000))
legend('V_n/V_A','V_n/c_{fast}','V_n/c_{slow}','1','Location','bestoutside','FontSize',16)
datetick('x','HH:MM')
ylim([0 2])
fshock.TileSpacing='none';

%% Plot dV & dVA
mu_0 = 4e-7*pi;
n0 = mean(np(subint)); VR0 = mean(vp(subint,1)); VT0 = mean(vp(subint,2)); VN0 = mean(vp(subint,3));
dBRp = BRp - mean(BRp); dBTp = BTp - mean(BTp); dBNp = BNp - mean(BNp);

dVAR = dBRp * 1e-9 / sqrt(mu_0*m_p*n0*1e6)/1e3; dVAT = dBTp * 1e-9 / sqrt(mu_0*m_p*n0*1e6)/1e3; dVAN = dBNp * 1e-9 / sqrt(mu_0*m_p*n0*1e6)/1e3;
dVR = vp(subint,1)-VR0; dVT = vp(subint,2)-VT0; dVN = vp(subint,3)-VN0;


f1=figure(1)
subplot(511)
thetaBR = acosd(BRp./Babsp);
plot(subepochp,dVAR)
title(date_str)
hold on
plot(subepochp,dVR)
plot([epochp(17097) epochp(17097)],[-50 50],'r--','LineWidth',1)
plot([epochp(17090) epochp(17090)],[-50 50],'r--','LineWidth',1)
plot([epochp(17110) epochp(17110)],[-50 50],'r--','LineWidth',1)
plot([epochp(17260) epochp(17260)],[-50 50],'r--','LineWidth',1)
plot([epochp(17280) epochp(17280)],[-50 50],'r--','LineWidth',1)
plot([epochp(17300) epochp(17300)],[-50 50],'r--','LineWidth',1)
plot([epochp(17400) epochp(17400)],[-50 50],'r--','LineWidth',1)
ylabel('dv_R (km/s)')
datetick('x','HH:MM')
subplot(512)
plot(subepochp,dVAT)
hold on
plot(subepochp,dVT)
ylabel('dv_T (km/s)')
datetick('x','HH:MM')
subplot(513)
plot(subepochp,dVAN)
hold on
plot(subepochp,dVN)
ylabel('dv_N (km/s)')
datetick('x','HH:MM')
subplot(514)
pcolor(tt(tt>=time_beg & tt <= time_end),f(fwave),ds_sigma_m(:,tt>=time_beg & tt <= time_end))
shading flat
set(gca,'yscale','log')
datetick('x','HH:MM')
ylabel('f (Hz)')
colormap jet
yyaxis right
plot(subepochp,thetaBR,'k')
ylabel('\theta_{BR}')
subplot(515)
plot(subepochp,np(subint))
ylabel('n_p (cm^{=3})')
yyaxis right
plot(subepochp,Tp(subint))
ylabel('T_p (eV)')
datetick('x','HH:MM')
%% load EAS data
easepoch = spdfcdfread(easpadefluxdir,'Variables',{'EPOCH'});
%easpadeflux = spdfcdfread(easpadefluxdir,'Variables',{'SWA_EAS_BM_Data'});
easpadeflux_arr = spdfcdfread(easpadefluxdir,'Variables',{'SWA_EAS2_Data'});
easEnergy_arr = spdfcdfread(easpadefluxdir,'Variables',{'SWA_EAS2_ENERGY'});
easEle = spdfcdfread(easpadefluxdir,'Variables',{'SWA_EAS_ELEVATION'});
easAzi = spdfcdfread(easpadefluxdir,'Variables',{'SWA_EAS_AZIMUTH'});
eas2RTN_arr = spdfcdfread(easpadefluxdir,'Variables',{'EAS2_TO_RTN'});
easdEle = spdfcdfread(easpadefluxdir,'Variables',{'SWA_EAS_ELEVATION_delta_upper'})+spdfcdfread(easpadefluxdir,'Variables',{'SWA_EAS_ELEVATION_delta_lower'});
easdAzi = spdfcdfread(easpadefluxdir,'Variables',{'SWA_EAS_AZIMUTH_delta_upper'})+spdfcdfread(easpadefluxdir,'Variables',{'SWA_EAS_AZIMUTH_delta_lower'});
%% calc eas pitch angle
easnum = numel(easepoch);
n_PAbin = 20;
strahlPA = zeros(easnum,64,n_PAbin);
PAnum = zeros(easnum,64,n_PAbin);
for etime = 1:easnum
easdEflux = easpadeflux_arr(:,:,:,etime); easEflux = easdEflux;
easEnergy = easEnergy_arr(etime,:);
%easEflux = easEflux./easEnergy';
eas2RTN = eas2RTN_arr(:,:,etime);
% for i = 1:64
%     for j = 1:32
%         for k = 1:16
%             easEflux(i,j,k)=easdEflux(i,j,k)*cosd(easEle(k))*easdEle(k)/180*pi*easdAzi(j)/180*pi;
%         end
%     end
% end
plotsub = epochmag >= easepoch(etime)-0*onesec & epochmag <= easepoch(etime)+100*onesec;
plotBR = mean(Bx(plotsub)); plotBT = mean(By(plotsub)); plotBN = mean(Bz(plotsub));
plotBeas = inv(eas2RTN)*[plotBR plotBT plotBN]'; plotBabs = sqrt(sum(plotBeas.^2));
[theta,phi]=meshgrid(easEle,easAzi);
PA = acosd((sind(theta).*cosd(phi)*plotBeas(1)+sind(theta).*sind(phi)*plotBeas(2)+cosd(theta)*plotBeas(3))/plotBabs);

for i = 1:64
    eastemp = squeeze(easEflux(i,:,:));
    eastemp(eastemp==0)=nan;
    for j = 1:n_PAbin
        strahlPA(etime,i,j)=sum(eastemp(PA>=(j-1)*(180/n_PAbin)&PA<j*(180/n_PAbin)),'omitnan');
        PAnum(etime,i,j)=numel(eastemp(PA>=(j-1)*(180/n_PAbin)&PA<j*(180/n_PAbin)));
    end
end
end
strahlPA(isnan(strahlPA))=0;
%% plot PA
subplot(611)
thetaBR = acosd(BRp./Babsp);
l1 = plot(subepochp,BRp);
hold on
l2 = plot(subepochp,sqrt(BRp.^2+BTp.^2+BNp.^2),'k');

title(date_str)
ylabel('B_R (nT)')
yyaxis right

plot(subepochp,vp(subint,1))
hold on 
%plot([epochp(17097) epochp(17097)],[300 400],'r--','LineWidth',1)
% plot([epochp(17090) epochp(17090)],[300 400],'r--','LineWidth',1)
% plot([epochp(17110) epochp(17110)],[300 400],'r--','LineWidth',1)
% plot([epochp(17260) epochp(17260)],[300 400],'r--','LineWidth',1)
% plot([epochp(17280) epochp(17280)],[300 400],'r--','LineWidth',1)
% plot([epochp(17300) epochp(17300)],[300 400],'r--','LineWidth',1)
legend([l1 l2],{'B_R','|B|'})
ylabel('v_R (km/s)')
datetick('x','HH:MM')
subplot(612)
plot(subepochp,BTp)
ylabel('B_T (nT)')
yyaxis right
plot(subepochp,vp(subint,2))
ylabel('v_T (km/s)')
datetick('x','HH:MM')
subplot(613)
plot(subepochp,BNp)
ylabel('B_N (nT)')
yyaxis right
plot(subepochp,vp(subint,3))
ylabel('v_N (km/s)')
datetick('x','HH:MM')
subplot(614)

pcolor(tt(tt>=time_beg & tt <= time_end),f(fwave),ds_sigma_m(:,tt>=time_beg & tt <= time_end))
title('\sigma_{mTN} (red = 1, blue = -1)')
shading flat
set(gca,'yscale','log')
datetick('x','HH:MM')
ylabel('f (Hz)')
yticks([0.01 0.1 0.5])
yticklabels({'0.01','0.1','1'})
colormap jet
yyaxis right
plot(subepochp,thetaBR,'k')
ylabel('\theta_{BR}')
ylim([0 180])
yticks([0 90 180])
yticklabels({'0','90','180'})
subplot(615)
plot(subepochp,wproton_arr(4,:))
ylabel('n_p (cm^{=3})')
yyaxis right
plot(subepochp,Tp(subint))
ylabel('T_p (eV)')
datetick('x','HH:MM')
subplot(616)
speedplot = speedbin <600 & speedbin > 200;
easint = easepoch>=time_beg & easepoch<=time_end;
pcolor(easepoch(easint),linspace(0,180,n_PAbin),log10(squeeze(sum(strahlPA(easint,20,:),2)))')
hold on
ylim([0 180])
%[maxnum,imaxmax]=max(eflux(subint,:)');
%plot(subepochp,speedbin(imaxmax),'w')
shading flat
title('>70 eV PA Eflux (eV m^-2 s^-1 eV^-1)')
ylabel('Pitch Angle')
datetick('x','HH:MM')

colormap jet
%% load vdf
epoch = spdfcdfread(cdfdir,'Variables',{'Epoch'});
energybin = double(spdfcdfread(cdfdir,'Variables',{'Energy'}));
speedbin = sqrt(energybin*eV*2/m_p)/1000; %km/s
azibin = double(spdfcdfread(cdfdir,'Variables',{'Azimuth'}));
elebin = double(spdfcdfread(cdfdir,'Variables',{'Elevation'}));
vdf = spdfcdfread(cdfdir,'Variables',{'vdf'});
vdfflag = spdfcdfread(cdfdir,'Variables',{'Info'});
%% plot vdf in RTN frame(ttime to select time)
savefigdir = 'E:\SolOquicklook\PAS_VDF\';
plotfitflag=1;
thetaBR = acosd(BRp./Babsp);
ttime = 17109; 
ttimebeg = 16020;dntime = 0;
Tcore = zeros(dntime+1,2); nforT=0;
for ttime = ttimebeg:ttimebeg+dntime
    nforT = nforT+1;
close all
pas2rtn_arr = spdfcdfread(cdfdir,'Variables',{'PAS_to_RTN'});
pas2rtn = pas2rtn_arr(:,:,ttime);
speedbin = sqrt(energybin*eV*2/m_p)/1000; %km/s
temp_vdf = double((vdf(:,:,:,ttime)));
[vsph,elesph,azisph] = meshgrid(speedbin,elebin,azibin);
vx = -cosd(elesph).*vsph.*cosd(azisph);
vy = -cosd(elesph).*vsph.*sind(azisph);
vz = -sind(elesph).*vsph;
vr = vx*pas2rtn(1,1)+vy*pas2rtn(1,2)+vz*pas2rtn(1,3);
vt = vx*pas2rtn(2,1)+vy*pas2rtn(2,2)+vz*pas2rtn(2,3);
vn = vx*pas2rtn(3,1)+vy*pas2rtn(3,2)+vz*pas2rtn(3,3);
maxVDF = max(temp_vdf,[],'all');
submax = find(temp_vdf==maxVDF);
f2 = figure;
subplot(2,3,1)

[vplot,phiplot]=meshgrid(speedbin,azibin/180*pi);
vxxyplot = vplot.*cos(phiplot);
vyxyplot = -vplot.*sin(phiplot)*pas2rtn(2,2);



size1d = [11*96*9 1];
finterp =  scatteredInterpolant(reshape(vr,size1d),reshape(vt,size1d),reshape(vn,size1d),reshape(temp_vdf,size1d));
finterp.ExtrapolationMethod =  'none';
ngrid = 50; vmax = 1000; dv = vmax/ngrid;
vxn = linspace(0,vmax,ngrid+1);
vyn = linspace(-vmax*sind(20),vmax*sind(40),ngrid+1);
vzn = linspace(-vmax/2,vmax/2,ngrid+1);
[vxnn,vynn,vznn]=meshgrid(vxn,vyn,vzn);
VDFxyz = finterp(vxnn,vynn,vznn);
VDFxyz(VDFxyz==0)=nan;

plotsub = epochmag >= epoch(ttime)-0*onesec & epochmag <= epoch(ttime)+1*onesec;
plotBR = mean(Bx(plotsub)); plotBT = mean(By(plotsub)); plotBN = mean(Bz(plotsub));
[maxnum,maxindex]=max(VDFxyz,[],'all','linear');
[maxi,maxj,maxk]= ind2sub(size(VDFxyz),maxindex);
contourf(vxn,vyn,log10(VDFxyz(:,:,maxk)),'ShowText','on')
f2.Position =  [281 160.2000 1.0856e+03 601.8000];
grid on
xlabel('V_R (km/s)')
ylabel('V_T (km/s)')
shading flat
colorbar
colormap jet
title(['log_{10} VDF (Cut of V_N =' num2str(vzn(maxk)) 'km/s)'])
hold on
plot(vxxyplot,vyxyplot,'k.')
plot(vp(ttime,1),vp(ttime,2),'wo','MarkerSize',10)
plot([vp(ttime,1) vp(ttime,1)+plotBR * 1e-9 / sqrt(mu_0*m_p*np(ttime)*1e6)/1e3],[vp(ttime,2) vp(ttime,2)+plotBT * 1e-9 / sqrt(mu_0*m_p*np(ttime)*1e6)/1e3],'w','LineWidth',2)
xlim([200 700])
ylim([-250 250])
axis square
subplot(2,3,2)

loglog(speedbin,eflux(ttime,:)'./energybin*m_p/eV*10,'x-')
xlabel('V (km/s)')
ylabel('Integrated 1D ion VDF (cm^{-3}km^{-1}s)')
title(datestr(epoch(ttime),'yyyymmdd HHMMSS') )
axis square
subplot(2,3,3)
plot(subepochp,thetaBR,'k')
ylabel('\theta_{BR}')
datetick('x','HH:MM')
ylim([0 180])
hold on
plot([epoch(ttime) epoch(ttime)],[0 180],'r--')
yyaxis right
plot(subepochp,np(subint))
ylabel('n_p (cm^{-3})')


% Transform to para perp frame
plotB0vec = [plotBR;plotBT;plotBN]; eb = plotB0vec/sqrt(sum(plotB0vec.^2));
k = cross([plotBR;plotBT;plotBN],[1;0;0]); k = k/sqrt(sum(k.^2));
k2 = cross(eb,k);
vbulkpara = double(vp(ttime,1)*eb(1)+vp(ttime,2)*eb(2)+vp(ttime,3)*eb(3));
vbulkperp1 = double(vp(ttime,1)*k(1)+vp(ttime,2)*k(2)+vp(ttime,3)*k(3));
vbulkperp2 = double(vp(ttime,1)*k2(1)+vp(ttime,2)*k2(2)+vp(ttime,3)*k2(3));
Trotate = [[eb(1) eb(2) eb(3)]; ...
    [k(1) k(2) k(3)];...
    [k2(1) k2(2) k2(3)]];
% xbphi = acos(eb(1));
% Trotate = [[cos(xbphi) -k(3)*sin(xbphi) k(2)*sin(xbphi)]; ...
%     [k(3)*sin(xbphi) k(2)^2+k(3)^2*cos(xbphi) k(2)*k(3)*(1-cos(xbphi))];...
%     [-k(2)*sin(xbphi) k(2)*k(3)*(1-cos(xbphi)) k(3)^2+k(2)^2*cos(xbphi)]];
vpara = Trotate(1,1)*vr+Trotate(1,2)*vt+Trotate(1,3)*vn-vbulkpara;
vperp1 = Trotate(2,1)*vr+Trotate(2,2)*vt+Trotate(2,3)*vn-vbulkperp1;
vperp2 = Trotate(3,1)*vr+Trotate(3,2)*vt+Trotate(3,3)*vn-vbulkperp2;
finterp2 =  scatteredInterpolant(reshape(vpara,size1d),reshape(vperp1,size1d),reshape(vperp2,size1d),reshape(temp_vdf,size1d));
finterp2.ExtrapolationMethod =  'none';
ngrid = 25; vmax = 500; dv = vmax/ngrid;
vxn = linspace(-vmax/2,vmax/2,ngrid+1);
vyn = linspace(-vmax/2,vmax/2,ngrid+1);
vzn = linspace(-vmax/2,vmax/2,ngrid+1);
[vxnn,vynn,vznn]=meshgrid(vxn,vyn,vzn);
VDFxyz = finterp2(vxnn,vynn,vznn);

[maxnum,maxindex]=max(VDFxyz,[],'all','linear');
[maxi,maxj,maxk]= ind2sub(size(VDFxyz),maxindex);
VDFxyz(VDFxyz==0)=nan;

size1d2 = [1 26*26*26];
VDF1D = reshape(VDFxyz,size1d2); vx1D = reshape(vxnn,size1d2); vy1D = reshape(vynn,size1d2); vz1D = reshape(vznn,size1d2);
vparaperp = [vx1D;vy1D;vz1D];
fitindex = log10(VDFxyz)>-8;
modelfun = @(b,v)(b(1)+(-(v(1,:)-b(2)).^2./(b(5).^2)-((v(2,:)-b(3)).^2+(v(3,:)-b(4)).^2)./(b(6).^2)));
beta0 = [-4,vxn(maxj),vyn(maxi),vzn(maxk),100,100];
beta = nlinfit(vparaperp(:,fitindex),log(VDF1D(fitindex)),modelfun,beta0);
vcore = (Trotate)\[beta(2)+vbulkpara;beta(3)+vbulkperp1;beta(4)+vbulkperp2];
Tcore(nforT,1)=beta(5); Tcore(nforT,2)=beta(6);
subplot(2,3,4)
contourf(vxn-beta(2),vyn-beta(3),log10(VDFxyz(:,:,maxk)),'ShowText','on')
hold on
theta_temp = linspace(0,2*pi,50); 
if plotfitflag == 1
A = 0.2; Tvx = A*(beta(5)*cos(theta_temp)); Tvy = A*(beta(6)*sin(theta_temp));
plot(Tvx,Tvy,'y--','LineWidth',1)
A = 0.4; Tvx = A*(beta(5)*cos(theta_temp)); Tvy = A*(beta(6)*sin(theta_temp));
plot(Tvx,Tvy,'y--','LineWidth',1)
A = 0.6; Tvx = A*(beta(5)*cos(theta_temp)); Tvy = A*(beta(6)*sin(theta_temp));
plot(Tvx,Tvy,'y--','LineWidth',1)
end
xlabel('V_{//} (km/s)')
ylabel('V_{\perp1} (km/s)')
shading flat
colorbar
colormap jet
title(['log_{10}VDF (Cut of V_{\perp2}=' num2str(vzn(maxk)) 'km/s)'])
axis square
grid on
subplot(2,3,5)
contourf(vxn-beta(2),vzn-beta(4),squeeze(log10(VDFxyz(maxi,:,:)))','ShowText','on')
hold on
theta_temp = linspace(0,2*pi,50); 
if plotfitflag == 1
A = 0.2; Tvx = A*(beta(5)*cos(theta_temp)); Tvy = A*(beta(6)*sin(theta_temp));
plot(Tvx,Tvy,'y--','LineWidth',1)
A = 0.4; Tvx = A*(beta(5)*cos(theta_temp)); Tvy = A*(beta(6)*sin(theta_temp));
plot(Tvx,Tvy,'y--','LineWidth',1)
A = 0.6; Tvx = A*(beta(5)*cos(theta_temp)); Tvy = A*(beta(6)*sin(theta_temp));
plot(Tvx,Tvy,'y--','LineWidth',1)
end
xlabel('V_{//} (km/s)')
ylabel('V_{\perp2} (km/s)')
shading flat
colorbar
colormap jet
title(['log_{10}VDF (Cut of V_{\perp1}=' num2str(vyn(maxi)) 'km/s)'])
axis square
grid on
subplot(2,3,6)
[maxnum,maxindex] = max(sum(sum(VDFxyz,3,'omitnan'),1));
contourf(vyn-beta(3),vzn-beta(4),squeeze(log10(VDFxyz(:,maxj,:)))','ShowText','on')
hold on
theta_temp = linspace(0,2*pi,50); 
if plotfitflag == 1
A = 0.2; Tvx = A*(beta(6)*cos(theta_temp)); Tvy = A*(beta(6)*sin(theta_temp));
plot(Tvx,Tvy,'y--','LineWidth',1)
A = 0.4; Tvx = A*(beta(6)*cos(theta_temp)); Tvy = A*(beta(6)*sin(theta_temp));
plot(Tvx,Tvy,'y--','LineWidth',1)
A = 0.6; Tvx = A*(beta(6)*cos(theta_temp)); Tvy = A*(beta(6)*sin(theta_temp));
plot(Tvx,Tvy,'y--','LineWidth',1)
end
xlabel('V_{\perp1} (km/s)')
ylabel('V_{\perp2} (km/s)')
shading flat
colorbar
colormap jet
title(['log_{10}VDF (Cut of V_{//}=' num2str(vxn(maxj)) 'km/s)'])
axis square
grid on
%close all
saveas(f2,[savefigdir datestr(epoch(ttime),'yyyymmdd_HHMMSS') '.png'])
end

f3 = figure(3);
plot(epoch(ttimebeg:ttimebeg+dntime),Tcore,'x','MarkerSize',20)
ylabel('v_{th} (km/s)')
ylim([0 200])
ttimesubindex = find(epoch==subepochp(1));
yyaxis right
plot(subepochp((ttimebeg:ttimebeg+dntime)-ttimesubindex+1),thetaBR((ttimebeg:ttimebeg+dntime)-ttimesubindex+1),'k')
ylabel('\theta_{BR}')
datetick('x','HH:MM:SS')
title([datestr(epoch(ttimebeg),'yyyymmdd HHMMSS') '\_' datestr(epoch(ttimebeg+dntime),'HHMMSS')])
legend('w_{//}','w_\perp','\theta_{BR}','Location','southeast')
saveas(f3,[savefigdir datestr(epoch(ttimebeg),'yyyymmddHHMMSS') '_' datestr(epoch(ttimebeg+dntime),'HHMMSS') '.png'])