clear;
fs = 20;
nhfile = 'Data/mar9_nh_antinu.mat'; % Data/feb2018_nh.mat
ihfile = 'Data/mar9_ih_antinu.mat'; % Data/mar2018_ih.mat
nhfile2 = 'Data/feb2018_nh.mat';
ihfile2 = 'Data/mar2018_ih.mat';
emuee_nh = load('Data/emuee_nh.mat');
emuee_ih = load('Data/emuee_ih.mat');
dune_nh = load(nhfile);     dune_nh2 = load(nhfile2);
dune_ih = load(ihfile);     dune_ih2 = load(ihfile2);
exp_nh = load(nhfile);      exp_nh2 = load(nhfile2);
exp_ih = load(ihfile);      exp_ih2 = load(ihfile2);
nonunit_nh = load(nhfile);  nonunit_nh2 = load(nhfile2);
nonunit_ih = load(ihfile);  nonunit_ih2 = load(ihfile2);
% dune_nh = load('Data/NH_data.mat');
% dune_ih = load('Data/IH_data.mat');
% exp_nh = load('Data/NH_data.mat');
% exp_ih = load('Data/IH_data.mat');
% nonunit_nh = load('Data/NH_data.mat');
% nonunit_ih = load('Data/IH_data.mat');
m = 0.0:0.002:0.20; % Lightest neutrino mass range
cflv_nh = emuee_nh.emu_ee;
cflv_ih = emuee_ih.emu_ee;
dune_nh = dune_nh.minValuesDune;    dune_nh2 = dune_nh2.minValuesDune;    
dune_ih = dune_ih.minValuesDune;    dune_ih2 = dune_ih2.minValuesDune;    
exp_nh = exp_nh.minValuesExp;       exp_nh2 = exp_nh2.minValuesExp;
exp_ih = exp_ih.minValuesExp;       exp_ih2 = exp_ih2.minValuesExp;
nonunit_nh = nonunit_nh.minValuesNonunit;
nonunit_ih = nonunit_ih.minValuesNonunit;
nonunit_nh2 = nonunit_nh2.minValuesNonunit;
nonunit_ih2 = nonunit_ih2.minValuesNonunit;
exp_nh = exp_nh-12; dune_nh = dune_nh-12; nonunit_nh = nonunit_nh-12;
exp_ih = exp_ih-12; dune_ih = dune_ih-12; nonunit_ih = nonunit_ih-12;
exp_nh2 = exp_nh2-12; dune_nh2 = dune_nh2-12; nonunit_nh2 = nonunit_nh2-12;
exp_ih2 = exp_ih2-12; dune_ih2 = dune_ih2-12; nonunit_ih2 = nonunit_ih2-12;

figure; % ANTI-NU
area(m,10.^(exp_nh),'FaceColor','y');
hold on;
area(m,10.^(dune_nh),'FaceColor', 'g'); % Grey = [0.5 0.5 0.5]
area(m,10.^(nonunit_nh),'FaceColor', 'c');
%area(m,cflv_nh,'FaceColor','r');
%cflv_nhplot = plot(m,cflv_nh,'-k');
%cflv_ihplot = plot(m,cflv_ih,'--k');
dune_nhplot = plot(m,10.^(dune_nh),'-k','LineWidth',2);
exp_nhplot = plot(m,10.^(exp_nh),'-k','LineWidth',2);
nonunit_nhplot = plot(m,10.^(nonunit_nh),'-k','LineWidth',2');
dune_ihplot = plot(m,10.^(dune_ih),'--k','LineWidth',2);
exp_ihplot = plot(m,10.^(exp_ih),'-.k','LineWidth',2);
nonunit_ihplot = plot(m,10.^(nonunit_ih),':k','LineWidth',2');
%plottables = [nonunit_nhplot nonunit_ihplot];
Beautify(false,0,1e13/1e12,fs);
text(0.07,7.5,'Excluded','FontSize',fs)
duneText = text(0.06,3.5,'DUNE coverage','FontSize',fs); set(duneText,'Rotation',-13);
text2 = text(0.06,1.42,'NSI not from nonunitarity','FontSize',fs); set(text2, 'Rotation', -3);
text3 = text(0.01,0.5,'Nonunitarity','FontSize',fs);
text(0.14,9,'anti-\nu','FontSize',fs+20);
%text3 = text(0.16,9,'NH','FontSize',fs+20);
%text(0.06,0.5,'Allowed by nonunitarity','FontSize',fs)
%text(0.07,9.3,'Allowed by CFLV','FontSize',fs)
set(gca,'XMinorTick','on','YMinorTick','on')
set(gca,'LineWidth',2,'TickLength',[0.025 0.025]);
grid on;
%legend(plottables,{'{\fontsize{15}NH}','{\fontsize{15}IH}'},'Location','NorthEast');

figure; % NU
area(m,10.^(exp_nh2),'FaceColor','y');
hold on;
area(m,10.^(dune_nh2),'FaceColor', 'g'); % Grey = [0.5 0.5 0.5]
area(m,10.^(nonunit_nh2),'FaceColor', 'c');
%area(m,cflv_nh,'FaceColor','r');
%cflv_nhplot = plot(m,cflv_nh,'-k');
%cflv_ihplot = plot(m,cflv_ih,'--k');
dune_nhplot = plot(m,10.^(dune_nh2),'-k','LineWidth',2);
exp_nhplot = plot(m,10.^(exp_nh2),'-k','LineWidth',2);
nonunit_nhplot = plot(m,10.^(nonunit_nh2),'-k','LineWidth',2');
dune_ihplot = plot(m,10.^(dune_ih2),'--k','LineWidth',2);
exp_ihplot = plot(m,10.^(exp_ih2),'-.k','LineWidth',2);
nonunit_ihplot = plot(m,10.^(nonunit_ih2),':k','LineWidth',2');
%plottables = [nonunit_nhplot nonunit_ihplot];
Beautify(false,0,1e13/1e12,fs);
text(0.07,7.5,'Excluded','FontSize',fs)
duneText = text(0.06,3.5,'DUNE coverage','FontSize',fs); set(duneText,'Rotation',-13);
text2 = text(0.06,1.42,'NSI not from nonunitarity','FontSize',fs); set(text2, 'Rotation', -3);
text3 = text(0.01,0.5,'Nonunitarity','FontSize',fs);
text(0.18,9,'\nu','FontSize',fs+20);
%text3 = text(0.16,9,'NH','FontSize',fs+20);
%text(0.06,0.5,'Allowed by nonunitarity','FontSize',fs)
%text(0.07,9.3,'Allowed by CFLV','FontSize',fs)
set(gca,'XMinorTick','on','YMinorTick','on')
set(gca,'LineWidth',2,'TickLength',[0.025 0.025]);
grid on;