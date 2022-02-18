function yy=ga_cba_RTN_fitness(c)
load tempdata.mat x11 x22 x33 yy eb eperp1 eperp2;
x1 = x11;
x2 = x22;
x3 = x33;
yt = yy;
% c(18) =
% [nc,nb,na,...
%  vcthp,vcthz,vbthp,vbthz,vathp,vathz,...
%  vcswpara,vcswperp1,vcswperp2,...
%  vbswpara,vbswperp1,vbswperp2,...
%  vcaswpara,vaswperp1,vaswperp2]
vpara = (x1-c(10)).*eb(1)+(x2-c(11)).*eb(2)+(x3-c(12)).*eb(3);
vperp1 = x1.*eperp1(1)+x2.*eperp1(2)+x3.*eperp1(3);
vperp2 = x1.*eperp2(1)+x2.*eperp2(2)+x3.*eperp2(3);
yf = 1.e-3*c(1)/(c(4)^2*c(5)*pi^1.5).*exp(-((x1-c(10)).*eb(1)+(x2-c(11)).*eb(2)+(x3-c(12)).*eb(3)).^2/c(5)^2).*exp(-((x1-c(10)).*eperp1(1)+(x2-c(11)).*eperp1(2)+(x3-c(12)).*eperp1(3)).^2/c(4)^2).*exp(-((x1-c(10)).*eperp2(1)+(x2-c(11)).*eperp2(2)+(x3-c(12)).*eperp2(3)).^2/c(4)^2)+...
     1.e-3*c(2)/(c(6)^2*c(7)*pi^1.5).*exp(-((x1-c(13)).*eb(1)+(x2-c(14)).*eb(2)+(x3-c(15)).*eb(3)).^2/c(7)^2).*exp(-((x1-c(13)).*eperp1(1)+(x2-c(14)).*eperp1(2)+(x3-c(15)).*eperp1(3)).^2/c(6)^2).*exp(-((x1-c(13)).*eperp2(1)+(x2-c(14)).*eperp2(2)+(x3-c(15)).*eperp2(3)).^2/c(6)^2)+...
     1.e-3*c(3)/(c(8)^2*c(9)*pi^1.5).*exp(-((x1-c(16)).*eb(1)+(x2-c(17)).*eb(2)+(x3-c(18)).*eb(3)).^2/c(9)^2).*exp(-((x1-c(16)).*eperp1(1)+(x2-c(17)).*eperp1(2)+(x3-c(18)).*eperp1(3)).^2/c(8)^2).*exp(-((x1-c(16)).*eperp2(1)+(x2-c(17)).*eperp2(2)+(x3-c(18)).*eperp2(3)).^2/c(8)^2);

%yy = sum((log10(yf)-log10(yt).^2/length(yt)));
end