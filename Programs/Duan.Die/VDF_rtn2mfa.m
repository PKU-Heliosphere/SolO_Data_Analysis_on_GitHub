function [VDFmfa,vxn,vxnn,vyn,vynn,vzn,vznn] = VDF_rtn2mfa(VDFrtn,vr,vt,vn,vrmax_vdf,vtmax_vdf,vnmax_vdf,eb,eperp1,eperp2)
%UNTITLED 此处显示有关此函数的摘要
is_vdf = VDFrtn>0;%&vr>520;
vdf0 = VDFrtn(is_vdf); vr0 = vr(is_vdf)-vrmax_vdf; vt0 = vt(is_vdf)-vtmax_vdf; vn0 = vn(is_vdf)-vnmax_vdf;
vpara0 = eb(1)*vr0+eb(2)*vt0+eb(3)*vn0;
vperp10= eperp1(1)*vr0+eperp1(2)*vt0+eperp1(3)*vn0;
vperp20 = eperp2(1)*vr0+eperp2(2)*vt0+eperp2(3)*vn0;
finterp = scatteredInterpolant(vpara0,vperp10,vperp20,vdf0);
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
%     vrnn = eb(1)*vxnn+eperp1(1)*vynn+eperp2(1)*vznn;
%     vtnn = eb(2)*vxnn+eperp1(2)*vynn+eperp2(2)*vznn;
%     vnnn = eb(3)*vxnn+eperp1(3)*vynn+eperp2(3)*vznn;
    VDFmfa = finterp(vxnn,vynn,vznn);  
end

