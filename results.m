clear
clc

% % Scordelis Lo-roof
% cps_2 = [4,6,10,18];
% cps_3 = [5,7,11,19];
% 
% % deg 2
% deg2 = [0.0715 0.2319 0.2969 0.3030];
% 
% % deg 3
% deg3 = [0.2972 0.3028 0.3034 0.3036];
% 
% % deg 2 ANS standard
% deg2_ANS = [0.2692 0.2994 0.3031 0.3034];
% 
% % deg 2 ANS new
% deg2_ANS_new = [0.3307 0.3249 0.3094 0.3052];
% 
% figure,
% plot(cps_2,deg2, '--x');
% hold on;
% plot(cps_3,deg3, '-o');
% hold on;
% plot(cps_2,deg2_ANS, '-v');
% hold on;
% plot(cps_2,deg2_ANS_new, '-o');
% hold on
% plot([0,20],[0.3024, 0.3024],'k', 'LineWidth',2);
% grid on;
% legend('p=2', 'p=3','p=2 ANS standard','p=2 ANS new','Location', 'northwest');
% xlabel('Control points per sides')
% ylabel('Vertical displacement')
% hold off;

% % Results Pinched
% cps_2 = [4,6,10,18,34];
% cps_3 = [5,7,11,19,35];
% 
% % deg 2
% deg2 = 1.0e-04*[0.0069   0.0164   0.0639   0.1446   0.1792];
% 
% % deg 3
% deg3 = 1.0e-04 *[0.0122   0.0534   0.1483   0.1815   0.1844];
% 
% % deg 2 ANS standard
% deg2_ANS = 1.0e-04*[0.0077   0.0230   0.0956   0.1714   0.1832];
% 
% % deg 2 ANS new
% deg2_ANS_new = 1.0e-04*[0.0160   0.0420   0.1334   0.1898   0.1902];
% 
% figure,
% plot(cps_2,deg2, '--x');
% hold on;
% plot(cps_3,deg3, '-o');
% hold on;
% plot(cps_2,deg2_ANS, '-v');
% hold on;
% plot(cps_2,deg2_ANS_new, '-o');
% hold on
% plot([0, 35],[1.8248e-5, 1.8248e-5],'k', 'LineWidth',2);
% grid on;
% legend('p=2', 'p=3','p=2 ANS standard','p=2 ANS new','Location', 'northwest');
% xlabel('Control points per sides')
% ylabel('Vertical displacement')
% hold off;


% Hemisphere
cps_2 = [4,6,10,18,34];
cps_3 = [5,7,11,19,35];

deg2 = [0.0001 0.0009 0.0131 0.0665 0.08989];
deg3 = [0.0013 0.0321 0.0879 0.0920 0.09245];
deg2_ANS = [0.0002 0.0078 0.0757 0.0909 0.0923];
deg2_ANS_new = [0.0005 0.0192 0.0858 0.0922 0.0927];

figure,
plot(cps_2,deg2, '--x');
hold on;
plot(cps_3,deg3, '-o');
hold on;
plot(cps_2,deg2_ANS, '-v');
hold on;
plot(cps_2,deg2_ANS_new, '-o');
hold on
plot([0, 35],[0.0924, 0.0924],'k', 'LineWidth',2);
grid on;
legend('p=2', 'p=3','p=2 ANS standard','p=2 ANS new','Location', 'northwest');
xlabel('Control points per sides')
ylabel('Radial displacement')
hold off;