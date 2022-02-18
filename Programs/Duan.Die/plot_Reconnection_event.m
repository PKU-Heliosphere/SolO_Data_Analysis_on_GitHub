Va_LMN = [BL;BM;BN]*walpha_arr(5:7,:)/sqrt(2);
Vc_LMN = [BL;BM;BN]*coeff_arr(:,10:12)';
figure(3)
subplot(2,1,1)
plot(epoch(sub_plot),walpha_arr(4,:))
yyaxis right
plot(epoch(sub_plot),coeff_arr(:,1))
legend('n_\alpha','n_c')
datetick('x','HH:MM')
xticklabels([])
subplot(2,1,2)
plot(epochmag_plot,BLMN(1,:))
datetick('x','HH:MM')
title('B_L (nT)')
figure(1)
subplot(4,1,1)
plot(epoch(sub_plot),Va_LMN(1,:))
hold on
plot(epoch(sub_plot),Vc_LMN(1,:))
datetick('x','HH:MM')
xticklabels([])
title('V_L (km/s)')


subplot(4,1,2)
plot(epoch(sub_plot),Va_LMN(2,:))

hold on
plot(epoch(sub_plot),Vc_LMN(2,:))
datetick('x','HH:MM')
xticklabels([])
title('V_M (km/s)')

subplot(4,1,3)
plot(epoch(sub_plot),Va_LMN(3,:))
hold on
plot(epoch(sub_plot),Vc_LMN(3,:))
datetick('x','HH:MM')
xticklabels([])
title('V_N (km/s)')

subplot(4,1,4)
plot(epochmag_plot,BLMN(1,:))
hold on
plot(epochmag_plot,vecnorm(BLMN,2,1))
datetick('x','HH:MM')
title('B (nT)')
%%
figure(2)
m_p = 1.673e-27;
eV = 1.602e-19;
nc = coeff_arr(:,1);
wc = coeff_arr(:,4:5);
Tc = wc.^2*2*m_p*1e6/eV;
subplot(4,1,1)
plot(epoch(sub_plot),Tc(:,2))
hold on
plot(epoch(sub_plot),Tc(:,1))
yyaxis right
plot(epoch(sub_plot),Vc_LMN(1,:),'k')
datetick('x','HH:MM')
xticklabels([])
title('T_c (eV)')
legend('T_{c//}','T_{c\perp}','V_L (km/s)')

subplot(4,1,2)
wa = walpha_arr(1:3,:)/sqrt(2);
Ta = wc.^2*2*4*m_p*1e6/eV;
plot(epoch(sub_plot),walpha_arr(1,:))
hold on
plot(epoch(sub_plot),walpha_arr(2,:))
plot(epoch(sub_plot),walpha_arr(3,:))
yyaxis right
plot(epoch(sub_plot),Va_LMN(1,:),'k')
datetick('x','HH:MM')
xticklabels([])
title('T_\alpha (eV)')
legend('T_{\alpha //}','T_{\alpha \perp 1}','T_{\alpha \perp 2}','V_{L\alpha}')
subplot(4,1,3)
plot(epoch(sub_plot),walpha_arr(8,:))
datetick('x','HH:MM')
title('Q_\alpha')
xticklabels([])

subplot(4,1,4)
plot(epochmag_plot,BLMN(1,:))
datetick('x','HH:MM')
title('B_L (nT)')