function [] = plotExperimentalData(EXPDATA)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

color11mM = [0.6350, 0.0780, 0.1840];           % Red
color_onlyliver = [0, 0.4470, 0.7410];          % Blue
          
marker='o';              

figureData=figure();
figureData.Position = [10 10 1300 800]; 

title ('Experimental data')

% Correct units
k_ins_nM=1/144; %Change insulin units to nM

subplot(2,2,1)
errorbar(EXPDATA.time{1}, EXPDATA.mean{1}, EXPDATA.SD{1},'Color', color11mM,...
    'LineStyle', 'n', 'Marker',marker,'MarkerFaceColor',color11mM, 'MarkerSize',8)
title('Glucose, liver-islet co-culture','FontSize',10)
xlabel('Time (h)');
ylabel('Glucose concentration (mM)','FontSize',20)
xticks([0 48 96 144 192 240 288 336])
set(gca,'TickDir','out','FontSize',20);
box off

subplot(2,2,2)
errorbar(EXPDATA.time{2}, EXPDATA.mean{2}*k_ins_nM, EXPDATA.SD{2}*k_ins_nM,'Color', color11mM,...
    'LineStyle', 'n', 'Marker',marker,'MarkerFaceColor',color11mM, 'MarkerSize',8)
title('Insulin, liver-islet co-culture','FontSize',10)
xlabel('Time (h)');
ylabel('Insulin concentration (nM)','FontSize',20)
xticks([0 48 96 144 192 240 288 336])
set(gca,'TickDir','out','FontSize',20);
box off

subplot(2,2,3)
errorbar(EXPDATA.time{3}, EXPDATA.mean{3}, EXPDATA.SD{3},'Color', color_onlyliver,...
    'LineStyle', 'n', 'Marker',marker,'MarkerFaceColor',color_onlyliver, 'MarkerSize',8)
title('Glucose, single-liver culture','FontSize',10)
xlabel('Time (h)');
ylabel('Glucose concentration (mM)','FontSize',20)
xticks([0 48 96 144 192 240 288 336])
set(gca,'TickDir','out','FontSize',20);
box off


end

