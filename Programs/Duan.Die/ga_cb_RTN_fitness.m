function yy=ga_cb_RTN_fitness(c)
load tempdata.mat x11 x22 x33 yy eb eperp1 eperp2;
x1 = x11;
x2 = x22;
x3 = x33;
yt = yy;
% c(18) =
% [nc,nb,...
%  vcthp,vcthz,vbthp,vbthz,...
%  vcswpara,vcswperp1,vcswperp2,...
%  vbswpara,vbswperp1,vbswperp2]
yf = 1.e-3*c(1)/(c(3)^2*c(4)*pi^1.5).*exp(-((x1-c(7)).*eb(1)+(x2-c(8)).*eb(2)+(x3-c(9)).*eb(3)).^2/c(4)^2).*exp(-((x1-c(7)).*eperp1(1)+(x2-c(8)).*eperp1(2)+(x3-c(9)).*eperp1(3)).^2/c(3)^2).*exp(-((x1-c(7)).*eperp2(1)+(x2-c(8)).*eperp2(2)+(x3-c(9)).*eperp2(3)).^2/c(3)^2)+...
     1.e-3*c(2)/(c(5)^2*c(6)*pi^1.5).*exp(-((x1-c(10)).*eb(1)+(x2-c(11)).*eb(2)+(x3-c(12)).*eb(3)).^2/c(6)^2).*exp(-((x1-c(10)).*eperp1(1)+(x2-c(11)).*eperp1(2)+(x3-c(12)).*eperp1(3)).^2/c(5)^2).*exp(-((x1-c(10)).*eperp2(1)+(x2-c(11)).*eperp2(2)+(x3-c(12)).*eperp2(3)).^2/c(5)^2);

yy = sum(sqrt((yf-yt).^2/length(yt)));
%yy = sum((log10(yf)-log10(yt)).^2/length(yt));
end