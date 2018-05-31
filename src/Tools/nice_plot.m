clear;clc;close all;

x0= 0.025; %m
v0= 0.000; %m/s
t = linspace(0,60,1000);
m = 1;     % kg
c = 0.1;   % N*sec/meter
k = 2.5;   % N/m

omega_n = sqrt(k/m)
zeta    = 0.5*c/m/omega_n

x = exp(-zeta*omega_n*t).*(x0*cos(sqrt(1-zeta^2)*omega_n*t) +...
   (v0+zeta*omega_n*x0)/sqrt(1-zeta^2)/omega_n*...
   (sin(sqrt(1-zeta^2)*omega_n*t)));

h=figure(1);
plot(t,x.*1000,'r-','LineWidth',4);
grid on;
FN = findall(h,'-property','FontName');
set(FN,'FontName','/usr/share/fonts/dejavu/DejaVuSerifCondensed.ttf');
FS = findall(h,'-property','FontSize');
set(FS,'FontSize',6);
title(["Response of an Underdamped Single Degree of\n"...
      "Freedom System Subjected to an Initial Excitation"],...
      'FontName','/usr/share/fonts/dejavu/DejaVuSerif-Italic.ttf',...
      'FontSize',8);
xlabel('Time (seconds)','FontSize',8);
ylabel('Displacement (mm)','FontSize',8);
L = legend(["x_0 = ", num2str(x0*1000,'%4.1f'), " mm    ",...
          "\nv_0 =  ", num2str(v0*1000,'%4.1f'), " mm/sec",...
          "\n  m =  ", num2str(m,'%4.1f'),       " kg    ",...
          "\n  c = ", num2str(c/1000,'%5.1e'),  " N-s/mm",...
          "\n  k = ", num2str(k/1000,'%5.1e'),  " N/mm  "])
FL1= findall(L,'-property','FontName');
FL2= findall(L,'-property','FontSize');
set(FL1,'FontName','/usr/share/fonts/msttcore/cour.ttf');
set(FL2,'FontSize',8);
H = 3; W = 4;
set(h,'PaperUnits','inches')
set(h,'PaperOrientation','portrait');
set(h,'PaperSize',[H,W])
set(h,'PaperPosition',[0,0,W,H])
print(h,'-dpng','-color','vib_plt4.png');