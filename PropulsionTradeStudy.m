%% Alex Nguyen
% Optimum Charateristics Trade Study
clc;clear all;close;
Rocket_EQN = @(MR) log(1/MR);
g = 9.81;
delv = 50;
mass_ratio = 0.886934;
x = [];
y = [];
for mass_ratio = .75:.005:.95
optimum_delv_c = Rocket_EQN(mass_ratio);
optimum_vc_c = fzero(@(x) (1/optimum_delv_c)*(exp(optimum_delv_c)-1)-1/2*x^2-1/2,1.01);
optimum_delv_vc = optimum_delv_c/optimum_vc_c;
optimum_vc = delv/optimum_delv_vc;

optimum_c = optimum_vc/optimum_vc_c;
optimum_isp = optimum_c/g;
x = [x,mass_ratio];
y = [y,optimum_isp];
end
% Propulsion System Parameters
plot(x,y,'k','linewidth',2)
hold on
grid on


name = 'MarCO';
mp_sys = 1.35679;%
m_sys = 12;%
Isp_sys = 40;%
tt = 30200;%
mdot_sys = mp_sys/tt; %
mpp_sys = 3.490;
plot((m_sys-mp_sys)/m_sys,Isp_sys,'b.','markersize',35)

name = 'MIPS';
mp_sys = .237003*2;%
m_sys = 12;%
Isp_sys = 40;%
tt = 9300;%
mdot_sys = mp_sys/tt; %
mpp_sys = .676*2;
plot((m_sys-mp_sys)/m_sys,Isp_sys,'g.','markersize',35)

name = 'VACCO Reaction  Control  Propulsion  Module';
P_sys = (9+12.6)/2*2;
mp_sys = .474006/(1.4)*2;%
m_sys = 12;%
F_sys = 20e-3*2;
Isp_sys = 40;%
mdot_sys = F_sys/Isp_sys/g*2; %
c_sys = Isp_sys/g*2;
mpp_sys = 1.244*2;
tt = mp_sys/mdot_sys;%
plot((m_sys-mp_sys)/m_sys,Isp_sys,'r.','markersize',35)
axis('tight')
% title('Optimal $I_{sp}$ for given payload mass fraction','interpreter','latex')
xlabel('Payload Mass Fraction')
legend('Optimal Line','MarCO','End Mounted MiPS','VACCO CPM','location','best')
ylabel('Optimal $I_{sp}$ [s]','interpreter','latex')
set(gca,'fontsize',18)
set(gcf,'paperorientation','landscape');
set(gcf,'paperunits','normalized');
set(gcf,'paperposition',[0 0 1 1]);
print(gcf,'-dpng','fig1.png');%output pdf name

fprintf('The optimal Isp for a mass payload ratio of %5.4f is %f\n',[mass_ratio,optimum_isp])
fprintf(['The ',name,' system has a Isp of %f\n'],Isp_sys)
fprintf(['The ',name,' system has a mass payload ratio of %f assuming a 12 kg system mass\n'],(m_sys-mp_sys)/m_sys)
fprintf('The expected total burn duration is 1600 s\n')
fprintf('The expected system total burn duration is %f [s]',tt)