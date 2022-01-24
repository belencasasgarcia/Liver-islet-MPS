function [] = plotExperimentalData(EXPDATA)


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

% Plot raw data to analyze systematic errors

data0h=[EXPDATA.G_raw(1,2:end) EXPDATA.G_raw(5,2:end) EXPDATA.G_raw(9,2:end) EXPDATA.G_raw(13,2:end)];
time0h=zeros(1,size(data0h,2));
data8h=[EXPDATA.G_raw(2,2:end) EXPDATA.G_raw(6,2:end) EXPDATA.G_raw(10,2:end) EXPDATA.G_raw(14,2:end)];
time8h=1*ones(1,size(data8h,2));
data24h=[EXPDATA.G_raw(3,2:end) EXPDATA.G_raw(7,2:end) EXPDATA.G_raw(11,2:end) EXPDATA.G_raw(15,2:end)];
time24h=2*ones(1,size(data24h,2));
data48h=[EXPDATA.G_raw(4,2:end) EXPDATA.G_raw(8,2:end) EXPDATA.G_raw(12,2:end) EXPDATA.G_raw(16,2:end)];
time48h=3*ones(1,size(data48h,2));

figure()

boxchart(time0h,data0h)
hold on
scatter(time0h, data0h,"filled",'jitter','on','JitterAmount',0.1,'MarkerFaceColor',[0, 0.4470, 0.7410]);
MarkerFaceAlpha = 0.5;
hold on
boxchart(time8h,data8h)
hold on
scatter(time8h, data8h,"filled",'jitter','on','JitterAmount',0.1,'MarkerFaceColor',[0.9290, 0.6940, 0.1250]);
MarkerFaceAlpha = 0.5;
hold on
boxchart(time24h,data24h)
hold on
scatter(time24h, data24h,"filled",'jitter','on','JitterAmount',0.1,'MarkerFaceColor',[0.4660, 0.6740, 0.1880]);
MarkerFaceAlpha = 0.5;
hold on
boxchart(time48h,data48h)
hold on
scatter(time48h, data48h,"filled",'jitter','on','JitterAmount',0.1,'MarkerFaceColor',[0.6350, 0.0780, 0.1840]);
MarkerFaceAlpha = 0.5;

xticks([0 1 2 3])
xticklabels({'0','8','24','48'})
xlabel('Time (h)');
ylabel('Glucose concentration (mM)','FontSize',20)
set(gca,'TickDir','out','FontSize',20);

end

