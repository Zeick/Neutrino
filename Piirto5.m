clear;
fs = 20;
emuee_nh = load('Data/emuee_nh.mat');
emuee_ih = load('Data/emuee_ih.mat');
dune_nh = load('Data/NH_data.mat');
dune_ih = load('Data/IH_data.mat');
exp_nh = load('Data/NH_data.mat');
exp_ih = load('Data/IH_data.mat');
nonunit_nh = load('Data/NH_data.mat');
nonunit_ih = load('Data/IH_data.mat');
m = 0.0:0.002:0.20; % Lightest neutrino mass range
cflv_nh = emuee_nh.emu_ee;
cflv_ih = emuee_ih.emu_ee;
dune_nh = dune_nh.minValuesDune;
dune_ih = dune_ih.minValuesDune;
exp_nh = exp_nh.minValuesExp;
exp_ih = exp_ih.minValuesExp;
nonunit_nh = nonunit_nh.minValuesNonunit;
nonunit_ih = nonunit_ih.minValuesNonunit;
area(m,10.^(exp_nh),'FaceColor','y');
hold on;
area(m,10.^(dune_nh),'FaceColor', 'g'); % Grey = [0.5 0.5 0.5]
area(m,10.^(nonunit_nh),'FaceColor', 'c');
%area(m,cflv_nh,'FaceColor','r');
%cflv_nhplot = plot(m,cflv_nh,'-k');
%cflv_ihplot = plot(m,cflv_ih,'--k');
dune_nhplot = plot(m,10.^(dune_nh),'-k');
dune_ihplot = plot(m,10.^(dune_ih),'--k');
exp_nhplot = plot(m,10.^exp_nh,'-k');
exp_ihplot = plot(m,10.^exp_ih,'--k');
nonunit_nhplot = plot(m,10.^(nonunit_nh),'-k');
nonunit_ihplot = plot(m,10.^(nonunit_ih),'--k');
plottables = [nonunit_nhplot nonunit_ihplot];
Beautify(false,0,1e13,fs);
%text(0.06,3.7e11,'Excluded','FontSize',fs)
%text(0.06,2.5e11,'NSI not from nonunitarity','FontSize',fs)
%text(0.06,0.5e11,'Allowed by nonunitarity','FontSize',fs)
%text(0.07,9.3,'Allowed by CFLV','FontSize',fs)

legend(plottables,{'{\fontsize{15}NH}','{\fontsize{15}IH}'},'Location','NorthEast');
