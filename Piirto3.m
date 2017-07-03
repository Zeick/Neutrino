clear;
fs = 20;
emuee_nh = load('Data/emuee_nh.mat');
emuee_ih = load('Data/emuee_ih.mat');
dune_nh = load('Data/DUNE_nh.mat');
dune_ih = load('Data/DUNE_ih.mat');
exp_nh = load('Data/Exp_NH.mat');
exp_ih = load('Data/Exp_IH.mat');
nonunit_nh = load('Data/Nonunit_nh.mat');
nonunit_ih = load('Data/Nonunit_ih.mat');
m = 0.0:0.002:0.20; % Lightest neutrino mass range
cflv_nh = emuee_nh.emu_ee;
cflv_ih = emuee_ih.emu_ee;
dune_nh = dune_nh.minValuesDune;
dune_ih = dune_ih.minValuesDune;
exp_nh = exp_nh.minValuesExp;
exp_ih = exp_ih.minValuesExp;
nonunit_nh = nonunit_nh.minValuesNonunit;
nonunit_ih = nonunit_ih.minValuesNonunit;
area(m,dune_nh,'FaceColor', 'g'); % Grey = [0.5 0.5 0.5]
hold on;
area(m,nonunit_nh,'FaceColor', 'c');
area(m,cflv_nh,'FaceColor','r');
cflv_nhplot = plot(m,cflv_nh,'-k');
cflv_ihplot = plot(m,cflv_ih,'--k');
dune_nhplot = plot(m,dune_nh,'-k');
dune_ihplot = plot(m,dune_ih,'--k');
%exp_nhplot = plot(m,exp_nh,'-r');
%exp_ihplot = plot(m,exp_ih,'--r');
nonunit_nhplot = plot(m,nonunit_nh,'-k');
nonunit_ihplot = plot(m,nonunit_ih,'--k');
plottables = [cflv_nhplot cflv_ihplot];
Beautify(false,9,12,fs);
text(0.07,11.7,'Excluded','FontSize',fs)
text(0.07,11.3,'NSI not from nonunitarity','FontSize',fs)
text(0.07,10.3,'Allowed by nonunitarity','FontSize',fs)
text(0.07,9.3,'Allowed by CFLV','FontSize',fs)

legend(plottables,{'{\fontsize{15}NH}','{\fontsize{15}IH}'},'Location','NorthEast');
