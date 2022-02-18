function [mom] = calc_3D_moment(f,y,x,z)
%CALC_3D_MOMENT 此处显示有关此函数的摘要
%   此处显示详细说明
mom2D = trapz(z,f,3);
mom1D = trapz(x,mom2D,2);
mom = trapz(y,mom1D);
end

