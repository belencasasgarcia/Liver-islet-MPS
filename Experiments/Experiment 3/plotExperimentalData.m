function [] = plotExperimentalData(EXPDATA)

color11mM = [0.6350, 0.0780, 0.1840];           % Red
color5p5mM = [0.4660, 0.6740, 0.1880];          % Blue
color2p8mM = [0.9290, 0.6940, 0.1250];	        % Yellow
          
marker='o';              

figureData=figure();
figureData.Position = [10 10 1300 800]; 

title ('Experimental data')

% Correct units
k_ins_nM=1/144; %Change insulin units to nM

subplot(3,2,1)
errorbar(EXPDATA.time{1}, EXPDATA.mean{1}, EXPDATA.SD{1},'Color', color11mM,...
    'LineStyle', 'n', 'Marker',marker,'MarkerFaceColor',color11mM, 'MarkerSize',8)
title('Glucose, liver-islet co-culture, hyperglycemia','FontSize',10)
xlabel('Time (h)');
ylabel('Glucose concentration (mM)','FontSize',20)
xticks([0 48 96 144 192 240 288 336])
set(gca,'TickDir','out','FontSize',20);
box off

subplot(3,2,2)
errorbar(EXPDATA.time{2}, EXPDATA.mean{2}*k_ins_nM, EXPDATA.SD{2}*k_ins_nM,'Color', color11mM,...
    'LineStyle', 'n', 'Marker',marker,'MarkerFaceColor',color11mM, 'MarkerSize',8)
title('Insulin, liver-islet co-culture,hyperglycemia', 'FontSize',10)
xlabel('Time (h)');
ylabel('Insulin concentration (nM)','FontSize',20)
xticks([0 48 96 144 192 240 288 336])
set(gca,'TickDir','out','FontSize',20);
box off

subplot(3,2,3)
errorbar(EXPDATA.time{5}, EXPDATA.mean{5}, EXPDATA.SD{5},'Color', color5p5mM,...
    'LineStyle', 'n', 'Marker',marker,'MarkerFaceColor',color5p5mM, 'MarkerSize',8)
title('Glucose, liver-islet co-culture, normoglycemia','FontSize',10)
xlabel('Time (h)');
ylabel('Glucose concentration (mM)','FontSize',20)
xticks([0 48 96 144 192 240 288 336])
set(gca,'TickDir','out','FontSize',20);
box off

subplot(3,2,4)
errorbar(EXPDATA.time{6}, EXPDATA.mean{6}*k_ins_nM, EXPDATA.SD{6}*k_ins_nM,'Color', color5p5mM,...
    'LineStyle', 'n', 'Marker',marker,'MarkerFaceColor',color5p5mM, 'MarkerSize',8)
title('Insulin, liver-islet co-culture, normoglycemia', 'FontSize',10)
xlabel('Time (h)');
ylabel('Insulin concentration (nM)','FontSize',20)
xticks([0 48 96 144 192 240 288 336])
set(gca,'TickDir','out','FontSize',20);
box off

subplot(3,2,5)
errorbar(EXPDATA.time{7}, EXPDATA.mean{7}, EXPDATA.SD{7},'Color', color2p8mM,...
    'LineStyle', 'n', 'Marker',marker,'MarkerFaceColor',color2p8mM, 'MarkerSize',8)
title('Glucose, liver-islet co-culture, hypoglycemia','FontSize',10)
xlabel('Time (h)');
ylabel('Glucose concentration (mM)','FontSize',20)
xticks([0 48 96 144 192 240 288 336])
set(gca,'TickDir','out','FontSize',20);
box off

subplot(3,2,6)
errorbar(EXPDATA.time{8}, EXPDATA.mean{8}*k_ins_nM, EXPDATA.SD{8}*k_ins_nM,'Color', color2p8mM,...
    'LineStyle', 'n', 'Marker',marker,'MarkerFaceColor',color2p8mM, 'MarkerSize',8)
title('Insulin, liver-islet co-culture, hypoglycemia', 'FontSize',10)
xlabel('Time (h)');
ylabel('Insulin concentration (nM)','FontSize',20)
xticks([0 48 96 144 192 240 288 336])
set(gca,'TickDir','out','FontSize',20);
box off

end

