%% This script calculate the uncertainty of the GA fitting method under the grids of SolO/PAS
m_p = 1.673e-27;%kg
eV = 1.602e-19;%J
k_b = 1.38e-23; 
%% load PAS grid
date_str = '20201014';
cdfdir = ['D:\SolOData\solo_l2_swa-pas-vdf_' date_str '_v01.cdf'];
epoch = spdfcdfread(cdfdir,'Variables',{'Epoch'});
energybin = double(spdfcdfread(cdfdir,'Variables',{'Energy'}));
azibin = double(spdfcdfread(cdfdir,'Variables',{'Azimuth'}));
elebin = double(spdfcdfread(cdfdir,'Variables',{'Elevation'}));
vdf = spdfcdfread(cdfdir,'Variables',{'vdf'});
pas2rtn_arr = spdfcdfread(cdfdir,'Variables',{'PAS_to_RTN'});

ttime = 29905;
epoch_tmp = epoch(ttime);
pas2rtn = pas2rtn_arr(:,:,ttime);
speedbin = sqrt(energybin*eV*2/m_p)/1000; %km/s
temp_vdf = double((vdf(:,:,:,ttime)));
[vsph,elesph,azisph] = meshgrid(speedbin,elebin,azibin);
vx = -cosd(elesph).*vsph.*cosd(azisph);
vy = -cosd(elesph).*vsph.*sind(azisph);
vz = sind(elesph).*vsph;
vrpas = vx*pas2rtn(1,1)+vy*pas2rtn(1,2)+vz*pas2rtn(1,3);
vtpas = vx*pas2rtn(2,1)+vy*pas2rtn(2,2)+vz*pas2rtn(2,3);
vnpas = vx*pas2rtn(3,1)+vy*pas2rtn(3,2)+vz*pas2rtn(3,3);

%%
ngridr = 41;  ngridt = 17;  ngridn = 17;
    vrmin = 300;  vrmax = 700; dvr = (vrmax-vrmin)/(ngridr-1);
    vtmin = -200;  vtmax = 200; dvt = (vtmax-vtmin)/(ngridt-1);% better to include vy = 0
    vnmin = -200;  vnmax = 200; dvn = (vnmax-vnmin)/(ngridn-1);% better to include vz = 0
    vrn = linspace(vrmin,vrmax,ngridr);
    vtn = linspace(vtmin,vtmax,ngridt);
    vnn = linspace(vnmin,vnmax,ngridn);
    [vrinterp,vtinterp,vninterp]=meshgrid(vrn,vtn,vnn);

%% calc VDF
%vr = vrpas; vt = vtpas; vn = vnpas;
vr = vrinterp; vt = vtinterp; vn = vninterp;

eb = [1,3,4];
eb = eb/norm(eb);
eperp2 = cross(eb,[0,1,0]);
eperp2 = eperp2/norm(eperp2);
eperp1 = cross(eperp2,eb);
eperp1 = eperp1/norm(eperp1);

nc = 13; vRc = 400; vTc =  20; vNc = -20; wparac = 50; wperpc  = 40;
nb = 1.6; vRb = 450; vTb = 0; vNb = 0; wparab = 90; wperpb = 90;
na = 1.4; vRa = 550; vTa = 0; vNa = 0; wparaa = 90; wperpa = 90;

VDFcore = 1.e-3*nc/(wperpc^2*wparac*pi^1.5).*exp(-((vr-vRc).*eb(1)+(vt-vTc).*eb(2)+(vn-vNc).*eb(3)).^2/wparac^2).*exp(-((vr-vRc).*eperp1(1)+(vt-vTc).*eperp1(2)+(vn-vNc).*eperp1(3)).^2/wperpc^2).*exp(-((vr-vRc).*eperp2(1)+(vt-vTc).*eperp2(2)+(vn-vNc).*eperp2(3)).^2/wperpc^2);
VDFbeam = 1.e-3*nb/(wperpb^2*wparab*pi^1.5).*exp(-((vr-vRb).*eb(1)+(vt-vTb).*eb(2)+(vn-vNb).*eb(3)).^2/wparab^2).*exp(-((vr-vRb).*eperp1(1)+(vt-vTb).*eperp1(2)+(vn-vNb).*eperp1(3)).^2/wperpb^2).*exp(-((vr-vRb).*eperp2(1)+(vt-vTb).*eperp2(2)+(vn-vNb).*eperp2(3)).^2/wperpb^2);
VDFalpha = 1.e-3*na/(wperpa^2*wparaa*pi^1.5).*exp(-((vr-vRa).*eb(1)+(vt-vTa).*eb(2)+(vn-vNa).*eb(3)).^2/wparaa^2).*exp(-((vr-vRa).*eperp1(1)+(vt-vTa).*eperp1(2)+(vn-vNa).*eperp1(3)).^2/wperpa^2).*exp(-((vr-vRa).*eperp2(1)+(vt-vTa).*eperp2(2)+(vn-vNa).*eperp2(3)).^2/wperpa^2);
%%
VDF_temp = VDFcore;
nc_mom = calc_3D_moment(VDF_temp,vtn,vrn,vnn);
vrc_mom = calc_3D_moment(VDF_temp.*vr,vtn,vrn,vnn)/nc_mom;
vtc_mom = calc_3D_moment(VDF_temp.*vt,vtn,vrn,vnn)/nc_mom;
vnc_mom = calc_3D_moment(VDF_temp.*vn,vtn,vrn,vnn)/nc_mom;

Pxx_mom = calc_3D_moment(VDF_temp.*(vr-vrc_mom).*(vr-vrc_mom),vtn,vrn,vnn);
wxx_mom = sqrt(Pxx_mom/nc_mom*2); Txx_mom = Pxx_mom*m_p/nc_mom*1e6/(eV);

Pyy_mom = calc_3D_moment(VDF_temp.*(vt-vtc_mom).*(vt-vtc_mom),vtn,vrn,vnn);
wyy_mom = sqrt(Pyy_mom/nc_mom*2); Tyy_mom = Pyy_mom*m_p/nc_mom*1e6/(eV);
Pzz_mom = calc_3D_moment(VDF_temp.*(vn-vnc_mom).*(vn-vnc_mom),vtn,vrn,vnn);
wzz_mom = sqrt(Pzz_mom/nc_mom*2); Tzz_mom = Pzz_mom*m_p/nc_mom*1e6/(eV);

vparaUpara = (vr-vrc_mom).*eb(1)+(vt-vtc_mom).*eb(2)+(vn-vnc_mom).*eb(3);
Ppara_mom = calc_3D_moment(VDF_temp.*vparaUpara.*vparaUpara,vtn,vrn,vnn);
wpara_mom = sqrt(Ppara_mom/nc_mom*2); Tpara_mom = Ppara_mom*m_p/nc_mom*1e6/(eV);

vperp1Uperp1 = (vr-vrc_mom).*eperp1(1)+(vt-vtc_mom).*eperp1(2)+(vn-vnc_mom).*eperp1(3);
Pperp1_mom = calc_3D_moment(VDF_temp.*vperp1Uperp1.*vperp1Uperp1,vtn,vrn,vnn);
wperp1_mom = sqrt(Pperp1_mom/nc_mom*2); Tperp1_mom = Pperp1_mom*m_p/nc_mom*1e6/(eV);

vperp2Uperp2 = (vr-vrc_mom).*eperp2(1)+(vt-vtc_mom).*eperp2(2)+(vn-vnc_mom).*eperp2(3);
Pperp2_mom = calc_3D_moment(VDF_temp.*vperp2Uperp2.*vperp2Uperp2,vtn,vrn,vnn);
wperp2_mom = sqrt(Pperp2_mom/nc_mom*2); Tperp2_mom = Pperp2_mom*m_p/nc_mom*1e6/(eV);

