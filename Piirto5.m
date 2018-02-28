clear;
fs = 20;
emuee_nh = load('Data/emuee_nh.mat');
emuee_ih = load('Data/emuee_ih.mat');
dune_nh = load('Data/feb2018_nh.mat');
dune_ih = load('Data/feb2018_ih.mat');
exp_nh = load('Data/feb2018_nh.mat');
exp_ih = load('Data/feb2018_ih.mat');
nonunit_nh = load('Data/feb2018_nh.mat');
nonunit_ih = load('Data/feb2018_ih.mat');
% dune_nh = load('Data/NH_data.mat');
% dune_ih = load('Data/IH_data.mat');
% exp_nh = load('Data/NH_data.mat');
% exp_ih = load('Data/IH_data.mat');
% nonunit_nh = load('Data/NH_data.mat');
% nonunit_ih = load('Data/IH_data.mat');
m = 0.0:0.002:0.20; % Lightest neutrino mass range
cflv_nh = emuee_nh.emu_ee;
cflv_ih = emuee_ih.emu_ee;
dune_nh = dune_nh.minValuesDune;
dune_ih = dune_ih.minValuesDune;
exp_nh = exp_nh.minValuesExp;
exp_ih = exp_ih.minValuesExp;
nonunit_nh = nonunit_nh.minValuesNonunit;
nonunit_ih = nonunit_ih.minValuesNonunit;
exp_nh = exp_nh-12; dune_nh = dune_nh-12; nonunit_nh = nonunit_nh-12;
exp_ih = exp_ih-12; dune_ih = dune_ih-12; nonunit_ih = nonunit_ih-12;

figure; % NH
area(m,10.^(exp_nh),'FaceColor','y');
hold on;
area(m,10.^(dune_nh),'FaceColor', 'g'); % Grey = [0.5 0.5 0.5]
area(m,10.^(nonunit_nh),'FaceColor', 'c');
%area(m,cflv_nh,'FaceColor','r');
%cflv_nhplot = plot(m,cflv_nh,'-k');
%cflv_ihplot = plot(m,cflv_ih,'--k');
dune_nhplot = plot(m,10.^(dune_nh),'-k','LineWidth',2);
exp_nhplot = plot(m,10.^exp_nh,'-k','LineWidth',2);
nonunit_nhplot = plot(m,10.^(nonunit_nh),'-k','LineWidth',2');
%plottables = [nonunit_nhplot nonunit_ihplot];
Beautify(false,0,1e13/1e12,fs);
text(0.06,7.5,'Excluded','FontSize',fs)
duneText = text(0.06,3.2,'DUNE coverage','FontSize',fs); set(duneText,'Rotation',-10);
text2 = text(0.06,1.9,'NSI not from nonunitarity','FontSize',fs); set(text2, 'Rotation', -5);
text3 = text(0.16,9,'NH','FontSize',fs+20);
%text(0.06,0.5,'Allowed by nonunitarity','FontSize',fs)
%text(0.07,9.3,'Allowed by CFLV','FontSize',fs)
set(gca,'XMinorTick','on','YMinorTick','on')
set(gca,'LineWidth',2,'TickLength',[0.025 0.025]);
grid on;
%legend(plottables,{'{\fontsize{15}NH}','{\fontsize{15}IH}'},'Location','NorthEast');

figure; % IH
area(m,10.^(exp_ih),'FaceColor','y');
hold on;
area(m,10.^(dune_ih),'FaceColor', 'g'); % Grey = [0.5 0.5 0.5]
area(m,10.^(nonunit_ih),'FaceColor', 'c');
dune_ihplot = plot(m,10.^(dune_ih),'-k','LineWidth',2);
exp_ihplot = plot(m,10.^exp_ih,'-k','LineWidth',2);
nonunit_ihplot = plot(m,10.^(nonunit_ih),'-k','LineWidth',2);
Beautify(false,0,1e13/1e12,fs);
text(0.06,7.5,'Excluded','FontSize',fs)
text3 = text(0.16,9,'IH','FontSize',fs+20);
set(gca,'XMinorTick','on','YMinorTick','on')
set(gca,'LineWidth',2,'TickLength',[0.025 0.025]);
grid on;