function [wpara,wperp1,wperp2,nc_mom,vrc_mom,vtc_mom,vnc_mom,Q] = calc_wpara_wperp(VDF,vtn,vrn,vnn,eb,eperp1,eperp2,n_M,T_or_w)
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明
m_p = 1.673e-27*n_M;%kg
eV = 1.602e-19;%J
k_b = 1.38e-23; 
VDF_temp = VDF;
[vr,vt,vn]=meshgrid(vrn,vtn,vnn);
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

Pxy_mom = calc_3D_moment(VDF_temp.*(vr-vrc_mom).*(vt-vtc_mom),vtn,vrn,vnn);
Pxz_mom = calc_3D_moment(VDF_temp.*(vr-vrc_mom).*(vn-vnc_mom),vtn,vrn,vnn);
Pyz_mom = calc_3D_moment(VDF_temp.*(vt-vtc_mom).*(vn-vnc_mom),vtn,vrn,vnn);

vparaUpara = (vr-vrc_mom).*eb(1)+(vt-vtc_mom).*eb(2)+(vn-vnc_mom).*eb(3);
Ppara_mom = calc_3D_moment(VDF_temp.*vparaUpara.*vparaUpara,vtn,vrn,vnn);
wpara = sqrt(Ppara_mom/nc_mom*2); Tpara_mom = Ppara_mom*m_p/nc_mom*1e6/(eV);

vperp1Uperp1 = (vr-vrc_mom).*eperp1(1)+(vt-vtc_mom).*eperp1(2)+(vn-vnc_mom).*eperp1(3);
Pperp1_mom = calc_3D_moment(VDF_temp.*vperp1Uperp1.*vperp1Uperp1,vtn,vrn,vnn);
wperp1 = sqrt(Pperp1_mom/nc_mom*2); Tperp1_mom = Pperp1_mom*m_p/nc_mom*1e6/(eV);

vperp2Uperp2 = (vr-vrc_mom).*eperp2(1)+(vt-vtc_mom).*eperp2(2)+(vn-vnc_mom).*eperp2(3);
Pperp2_mom = calc_3D_moment(VDF_temp.*vperp2Uperp2.*vperp2Uperp2,vtn,vrn,vnn);
wperp2 = sqrt(Pperp2_mom/nc_mom*2); Tperp2_mom = Pperp2_mom*m_p/nc_mom*1e6/(eV);

P12_mom = calc_3D_moment(VDF_temp.*vparaUpara.*vperp1Uperp1,vtn,vrn,vnn);
P13_mom = calc_3D_moment(VDF_temp.*vparaUpara.*vperp2Uperp2,vtn,vrn,vnn);
P23_mom = calc_3D_moment(VDF_temp.*vperp1Uperp1.*vperp2Uperp2,vtn,vrn,vnn);
Q = (P12_mom^2+P13_mom^2+P23_mom^2)/(Ppara_mom*Pperp1_mom+Ppara_mom*Pperp2_mom+Pperp1_mom*Pperp2_mom);

if T_or_w == 'T'
    wpara = Tpara_mom;
    wperp1 = Tperp1_mom;
    wperp2 = Tperp2_mom;
end
%I1 = Pxx_mom + Pyy_mom + Pzz_mom;
%I2 = Pxx_mom*Pyy_mom

end

