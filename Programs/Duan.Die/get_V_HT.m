function V_HT = get_V_HT(v,B)
%Get velocity of deHoffmann-Tellor frame
%   input: v(n*3), time series of velocity vectors
%          B(n*3), time series of magnetic field vectors
%   output: v_HT, velocity of deHoffmann-Tellor frame
    norm_B = sqrt(sum(B.^2,2));
    K = zeros(3,3);
    K_dot_v = zeros(3,3);
    for i = 1:3
        for j= 1:3
           K(i,j) = mean(norm_B.^2.*kroneckerDelta(sym(i),sym(j))-B(:,i).*B(:,j));
           K_dot_v(i,j) = mean((norm_B.^2.*kroneckerDelta(sym(i),sym(j))-B(:,i).*B(:,j)).*v(:,j));
        end
    end
    K_dot_v = sum(K_dot_v,2);
    V_HT = K^(-1)*K_dot_v;
    V_HT = [double(V_HT(1)), double(V_HT(2)), double(V_HT(3))];
end