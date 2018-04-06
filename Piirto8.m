clear;
fs = 20;
nhfile = 'Data/mar11_nh_antinu.mat';
ihfile = 'Data/mar11_ih_antinu.mat';
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

j=101; % 1 -> 0 eV, 51 -> 0.1 eV, 101 -> 0.2 eV 
figure;
m1 = (j-1)*0.002;
dune1_lin = 10.^dune_nh(j);            dune2_lin = 10.^dune_ih2(j);
exp1_lin = 10.^exp_nh(j);              exp2_lin = 10.^exp_ih2(j);
nonunit1_lin = 10.^nonunit_nh(j);      nonunit2_lin = 10.^nonunit_ih2(j);
mDelta = 0:0.1:10; % TeV                
dune1_lin = mDelta/dune1_lin;          dune2_lin = mDelta/dune2_lin;
exp1_lin = mDelta/exp1_lin;            exp2_lin = mDelta/exp2_lin;
nonunit1_lin = mDelta/nonunit1_lin;    nonunit2_lin = mDelta/nonunit2_lin;
hold on;
ylim([0 3]);
set(gca,'Color','c','FontSize',fs); % Background color of the plot
set(gcf,'color','w'); % Background color of the plot window
area(mDelta,nonunit2_lin,'FaceColor','g');
area(mDelta,dune2_lin,'FaceColor', 'y');
area(mDelta,exp2_lin,'FaceColor', 'w'); % Grey = [0.5 0.5 0.5]
%plot(mDelta, dune1_lin,'-k');
plot(mDelta, dune2_lin,'-k');
%plot(mDelta, exp1_lin,'-k');
plot(mDelta, exp2_lin,'-k');
%plot(mDelta, nonunit1_lin,'-k');
plot(mDelta, nonunit2_lin,'-k');
xlabel('M_{\Delta} (TeV)','FontSize',fs);
ylabel('\lambda_{\phi} (eV)','FontSize',fs);
text(8,0.5,'IH','FontSize',fs+20);
%title(['m_1 = ' num2str(m1) ' eV'],'FontSize',fs);
%legend('{\fontsize{15}Experimental limit}', '{\fontsize{15}DUNE coverage}','{\fontsize{15}Non-unitary limit}','Location','NorthEast');
hold off;
if j < 10
    filename = sprintf('Fig3_nh0%d.png',j);
else
    filename = sprintf('Fig3_nh%d.png',j);
end
%    saveas(gcf,filename);