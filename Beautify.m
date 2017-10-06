% (C) Timo Karkkainen 2017
% Helper-function. low = lowest point of y-axis visible, high = highest.
% Makes the plot more human-readable.
function Beautify(showLegend,low, high, fontsize)
set(gca,'FontSize',fontsize);
xlabel('m_1 (eV)','FontSize',fontsize);
%ylabel('log_{10}(M_{\Delta}/|\lambda_{\phi}|)','FontSize',fontsize);
ylabel('M_{\Delta}/|\lambda_{\phi}|  (\times10^{12})','FontSize',fontsize);
%title('M_{\Delta}/|\lambda_{\phi}| as a function of the lightest neutrino','FontSize',fontsize);
ylim([low high]);
if(showLegend)
    legend('{\fontsize{15}ee-\mu\mu}','{\fontsize{15}e\mu}','{\fontsize{15}e\tau}','{\fontsize{15}\mu\tau}','{\fontsize{15}\tau\tau -\mu\mu}','Location','NorthEast');
end