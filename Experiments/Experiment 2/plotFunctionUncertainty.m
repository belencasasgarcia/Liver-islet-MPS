function plotFunctionUncertainty(param, model, parIndex, param_opt, N_exchanges,statesIndex)

time_end=335.99;

loadExperimentalData;

COSTOPTIONS = [];
COSTOPTIONS.maxnumsteps = 1e8;
COSTOPTIONS.abstol = 1e-10;
COSTOPTIONS.reltol = 1e-10;
COSTOPTIONS.minstep = 0;
COSTOPTIONS.maxstep = inf;

%% Plotting settings

color11mM = [0.6350, 0.0780, 0.1840];           % Red
color_onlyliver = [0, 0.4470, 0.7410];          % Blue
color_onlyislets = [0.8500, 0.3250, 0.0980];	% Orange
          
marker='o';              

k_ins_nM=1/144;

modelname = char(model);

%simTime=repmat([0:0.01:(L_GTT-0.01)],N_exchanges,1)+[0:L_GTT:L_GTT*(N_exchanges-1)]';
%simTime=[0:0.01:time_end];

% Define simulation based on the number of media exchanges
simTime=[0:0.01:(L_GTT*N_exchanges-0.01)];
nTime=size(simTime,2);
time_intervals=boolean(zeros(N_exchanges,size(simTime,2)));

% Select time points corresponding to each media exchange "k"

time_intervals(1,1:size([0:0.01:47.99],2))=1;

for k=2:size(time_intervals,1)
    time_intervals(k,:)=circshift(time_intervals(k-1,:),size([0:0.01:47.99],2));
end


%% Simulate for liver-islet co-culture, hyperglycemia (11 mM)

% Account for offset in glucose and insulin measurements due to media
% exchanges

param(:,parIndex.i_delta_G_1)=param(:,parIndex.i_delta_G_1_11mM);
param(:,parIndex.i_delta_G_7)=0;
param(:,parIndex.i_delta_G_13)=param(:,parIndex.i_delta_G_13_11mM);

param(:,parIndex.i_delta_I_1)=param(:,parIndex.i_delta_I_1_11mM);
param(:,parIndex.i_delta_I_7)=0;
param(:,parIndex.i_delta_I_13)=param(:,parIndex.i_delta_I_13_11mM);

%simTime=reshape(simTime,1,size(simTime,2)*N_exchanges);

plotData_11mM = simulateCycles(param,parIndex,modelname,time_intervals,COSTOPTIONS,N_exchanges,k_ins_nM);

% Calculate maximal and minimal values of the simulations to represent uncertainties

plotData_11mM.Gmeasured_max(1:nTime)=max(plotData_11mM.Gmeasured,[],1);
plotData_11mM.Gmeasured_min(1:nTime)=min(plotData_11mM.Gmeasured,[],1);
plotData_11mM.Imeasured_max(1:nTime)=max(plotData_11mM.Imeasured,[],1);
plotData_11mM.Imeasured_min(1:nTime)=min(plotData_11mM.Imeasured,[],1);
plotData_11mM.Gliver_max(1:nTime)=max(plotData_11mM.Gliver,[],1);
plotData_11mM.Gliver_min(1:nTime)=min(plotData_11mM.Gliver,[],1);
plotData_11mM.Gislets_max(1:nTime)=max(plotData_11mM.Gislets,[],1);
plotData_11mM.Gislets_min(1:nTime)=min(plotData_11mM.Gislets,[],1);
plotData_11mM.Iliver_max(1:nTime)=max(plotData_11mM.Iliver,[],1);
plotData_11mM.Iliver_min(1:nTime)=min(plotData_11mM.Iliver,[],1);
plotData_11mM.Iislets_max(1:nTime)=max(plotData_11mM.Iislets,[],1);
plotData_11mM.Iislets_min(1:nTime)=min(plotData_11mM.Iislets,[],1);
plotData_11mM.U_ii_min(1:nTime)=min(plotData_11mM.U_ii,[],1);
plotData_11mM.U_ii_max(1:nTime)=max(plotData_11mM.U_ii,[],1);
plotData_11mM.U_id_min(1:nTime)=min(plotData_11mM.U_id,[],1);
plotData_11mM.U_id_max(1:nTime)=max(plotData_11mM.U_id,[],1);
plotData_11mM.Glucose_int_max(1:nTime) = max(plotData_11mM.Glucose_int,[],1);
plotData_11mM.Glucose_int_min(1:nTime) = min(plotData_11mM.Glucose_int,[],1);
plotData_11mM.SI_max(1:nTime) = max(plotData_11mM.SI,[],1);
plotData_11mM.SI_min(1:nTime) = min(plotData_11mM.SI,[],1);
plotData_11mM.Gslow_max(1:nTime) = max(plotData_11mM.Gslow,[],1);
plotData_11mM.Gslow_min(1:nTime) = min(plotData_11mM.Gslow,[],1);
plotData_11mM.Vislets_max(1:nTime) = max(plotData_11mM.Vislets,[],1);
plotData_11mM.Vislets_min(1:nTime) = min(plotData_11mM.Vislets,[],1);
plotData_11mM.Secretion_max(1:nTime) = max(plotData_11mM.Secretion,[],1);
plotData_11mM.Secretion_min(1:nTime) = min(plotData_11mM.Secretion,[],1);

% Simulate for the optimal parameter values

param_opt(parIndex.i_delta_G_1)=param_opt(parIndex.i_delta_G_1_11mM);
param_opt(parIndex.i_delta_G_7)=param_opt(parIndex.i_delta_G_7_11mM);
param_opt(parIndex.i_delta_G_13)=param_opt(parIndex.i_delta_G_13_11mM);

param_opt(parIndex.i_delta_I_1)=param_opt(parIndex.i_delta_I_1_11mM);
param_opt(parIndex.i_delta_I_7)=param_opt(parIndex.i_delta_I_7_11mM);
param_opt(parIndex.i_delta_I_13)=param_opt(parIndex.i_delta_I_13_11mM);

plotData_11mM_opt = simulateCycles(param_opt,parIndex,modelname,time_intervals,COSTOPTIONS,N_exchanges,k_ins_nM);

plotData_11mM.Gmeasured_opt = plotData_11mM_opt.Gmeasured;
plotData_11mM.Imeasured_opt = plotData_11mM_opt.Imeasured;
plotData_11mM.Gliver_opt = plotData_11mM_opt.Gliver;
plotData_11mM.Gislets_opt = plotData_11mM_opt.Gislets;
plotData_11mM.Iliver_opt = plotData_11mM_opt.Iliver;
plotData_11mM.Iislets_opt = plotData_11mM_opt.Iislets;
plotData_11mM.U_ii_opt = plotData_11mM_opt.U_ii;
plotData_11mM.U_id_opt = plotData_11mM_opt.U_id;

plotData_11mM.Glucose_int_opt = plotData_11mM_opt.Glucose_int;
plotData_11mM.SI_opt = plotData_11mM_opt.SI;
plotData_11mM.Gslow_opt = plotData_11mM_opt.Gslow;
plotData_11mM.Vislets_opt = plotData_11mM_opt.Vislets;
plotData_11mM.Secretion_opt = plotData_11mM_opt.Secretion;

% Disease progression variables

figure1=figure();
figure1.Position = [10 10 1500 1000]; 

% Time integral of excess glucose

subplot(2,3,1)

hold on

area1=[simTime,fliplr(simTime)];
area2=[plotData_11mM.Glucose_int_max,fliplr(plotData_11mM.Glucose_int_min)];
f=fill(area1,area2,color11mM,'LineStyle','none');
alpha(f,.2)

hold on
plot(simTime,plotData_11mM.Glucose_int_opt,'Color', color11mM,'Linewidth',2);
hold off

xlim([0 simTime(end)])
ylim([0 1200])
xticks([0 48 96 144 192 240 288 336])
xlabel('Time (h)','FontSize',20);
ylabel('AUC excess glucose (mM\cdoth)','FontSize',20)
set(gca,'TickDir','out','FontSize',20);
box off

title('Integral of excess glucose over time')

% Hepatic insulin sensitivity

subplot(2,3,2)

hold on

area1=[simTime,fliplr(simTime)];
area2=[plotData_11mM.SI_max,fliplr(plotData_11mM.SI_min)];
f=fill(area1,area2,color11mM,'LineStyle','none');
alpha(f,.2)

hold on
plot(simTime,plotData_11mM.SI_opt,'Color', color11mM,'Linewidth',2);
hold off

xlim([0 simTime(end)])
ylim([0 0.5])
xlabel('Time (h)','FontSize',20);
ylabel('Hepatic insulin resistance ^{L}/_{mIU\cdoth}','FontSize',20)
set(gca,'TickDir','out','FontSize',20);
box off

title('Hepatic insulin sensitivity')

% Daily average glucose

subplot(2,3,3)

hold on

area1=[simTime,fliplr(simTime)];
area2=[plotData_11mM.Gslow_max,fliplr(plotData_11mM.Gslow_min)];
f=fill(area1,area2,color11mM,'LineStyle','none');
alpha(f,.2)

hold on
plot(simTime,plotData_11mM.Gslow_opt,'Color', color11mM,'Linewidth',2);
hold off

xlim([0 simTime(end)])
ylim([5.5 7])
xticks([0 48 96 144 192 240 288 336])
xlabel('Time (h)','FontSize',20);
ylabel('Daily average glucose (mM)','FontSize',20)
set(gca,'TickDir','out','FontSize',20);
box off

title('Daily average glucose concentration')

%Beta cell volume

subplot(2,3,4)

hold on

area1=[simTime,fliplr(simTime)];
area2=[plotData_11mM.Vislets_max*1e9,fliplr(plotData_11mM.Vislets_min*1e9)];
f=fill(area1,area2,color11mM,'LineStyle','none');
alpha(f,.2)

hold on
plot(simTime,plotData_11mM.Vislets_opt*1e9,'Color', color11mM,'Linewidth',2);
hold off

xlim([0 simTime(end)])
ylim([8 20])
xticks([0 48 96 144 192 240 288 336])
xlabel('Time (h)','FontSize',20);
ylabel('\beta cell volume (nL)','FontSize',20)
set(gca,'TickDir','out','FontSize',20);
box off

title('Changes in \beta-cell volume')

%Beta cell insulin secretion capacity

subplot(2,3,5)

hold on

area1=[simTime,fliplr(simTime)];
area2=[plotData_11mM.Secretion_max*100,fliplr(plotData_11mM.Secretion_min*100)];
f=fill(area1,area2,color11mM,'LineStyle','none');
alpha(f,.2)

hold on
plot(simTime,plotData_11mM.Secretion_opt*100,'Color', color11mM,'Linewidth',2);
hold off

xlim([0 simTime(end)])
ylim([0 100])
xticks([0 48 96 144 192 240 288 336])
xlabel('Time (h)','FontSize',20);
ylabel('Insulin secretion capacity (% per max)','FontSize',20)
set(gca,'TickDir','out','FontSize',20);
box off

title('Changes in \beta-cell insulin secretion')

% Glucose and insulin in all compartments

figure2=figure();
figure2.Position = [10 10 1500 1000]; 

% Glucose in liver compartment

subplot(2,3,1)

hold on
area1=[simTime,fliplr(simTime)];
area2=[plotData_11mM.Gliver_max,fliplr(plotData_11mM.Gliver_min)];
f=fill(area1,area2,color_onlyliver,'LineStyle','none');
alpha(f,.4)

hold on
errorbar(EXPDATA.time{4}, EXPDATA.mean{4}, EXPDATA.SD{4},'Color', color_onlyliver,...
    'LineStyle', 'n', 'Marker',marker,'MarkerFaceColor',color_onlyliver, 'MarkerSize',8)

hold off
title('Liver compartment, liver-islet co-culture','FontSize',10)
xlabel('Time (h)');
ylabel('Glucose concentration (mM)','FontSize',20)
xlim([0 simTime(end)])
ylim([0 12])
xticks([0 48 96 144 192 240 288 336])
set(gca,'TickDir','out','FontSize',20);
box off

% Glucose in the islets compartment

subplot(2,3,2)

hold on
area1=[simTime,fliplr(simTime)];
area2=[plotData_11mM.Gislets_max,fliplr(plotData_11mM.Gislets_min)];
f=fill(area1,area2,color_onlyislets,'LineStyle','none');
alpha(f,.4)

hold on
errorbar(EXPDATA.time{5}, EXPDATA.mean{5}, EXPDATA.SD{5},'Color', color_onlyislets,...
    'LineStyle', 'n', 'Marker',marker,'MarkerFaceColor',color_onlyislets, 'MarkerSize',8)

hold off
title('Pancreas compartment, liver-islet co-culture','FontSize',10)
xlabel('Time (h)');
ylabel('Glucose concentration (mM)','FontSize',20)
xlim([0 simTime(end)])
ylim([0 12])
xticks([0 48 96 144 192 240 288 336])
set(gca,'TickDir','out','FontSize',20);
box off

% Measured glucose

subplot(2,3,3)

hold on

area1=[simTime,fliplr(simTime)];
area2=[plotData_11mM.Gmeasured_max,fliplr(plotData_11mM.Gmeasured_min)];
f=fill(area1,area2,color11mM,'LineStyle','none');
alpha(f,.4)

hold on
errorbar(EXPDATA.time{1}, EXPDATA.mean{1}, EXPDATA.SD{1},'Color', color11mM,...
    'LineStyle', 'n', 'Marker',marker,'MarkerFaceColor',color11mM, 'MarkerSize',8)
hold on
plot(simTime,plotData_11mM.Gmeasured_opt,'Color', color11mM,'Linewidth',2);
hold off

xlim([0 simTime(end)])
ylim([0 12])
xticks([0 48 96 144 192 240 288 336])
title('Measured, liver-islet co-culture','FontSize',10)
xlabel('Time (h)','FontSize',20);
ylabel('Glucose concentration (mM)','FontSize',20)
set(gca,'TickDir','out','FontSize',20);
box off

% Insulin in the liver compartment

subplot(2,3,4)

hold on
area1=[simTime,fliplr(simTime)];
area2=[plotData_11mM.Iliver_max,fliplr(plotData_11mM.Iliver_min)];
f=fill(area1,area2,color_onlyliver,'LineStyle','none');
alpha(f,.4)

hold on
errorbar(EXPDATA.time{6}, EXPDATA.mean{6}*k_ins_nM, EXPDATA.SD{6}*k_ins_nM,'Color', color_onlyliver,...
    'LineStyle', 'n', 'Marker',marker,'MarkerFaceColor',color_onlyliver, 'MarkerSize',8)

hold off
title('Liver compartment, liver-islet co-culture','FontSize',10)
xlabel('Time (h)');
ylabel('Insulin concentration (nM)','FontSize',20)
xlim([0 simTime(end)])
xticks([0 48 96 144 192 240 288 336])
set(gca,'TickDir','out','FontSize',20);
box off

% Insulin in the islets compartment

subplot(2,3,5)

hold on
area1=[simTime,fliplr(simTime)];
area2=[plotData_11mM.Iislets_max,fliplr(plotData_11mM.Iislets_min)];
f=fill(area1,area2,color_onlyislets,'LineStyle','none');
alpha(f,.4)

hold on
errorbar(EXPDATA.time{7}, EXPDATA.mean{7}*k_ins_nM, EXPDATA.SD{7}*k_ins_nM,'Color', color_onlyislets,...
    'LineStyle', 'n', 'Marker',marker,'MarkerFaceColor',color_onlyislets, 'MarkerSize',8)

hold off
title('Pancreas compartment, liver-islet co-culture','FontSize',10)
xticks([0 48 96 144 192 240 288 336])
xlabel('Time (h)');
ylabel('Insulin concentration (nM)','FontSize',20)
xlim([0 simTime(end)])
set(gca,'TickDir','out','FontSize',20);
box off

% Measured insulin

subplot(2,3,6)

hold on

area1=[simTime,fliplr(simTime)];
area2=[plotData_11mM.Imeasured_max,fliplr(plotData_11mM.Imeasured_min)];
f=fill(area1,area2,color11mM,'LineStyle','none');
alpha(f,.2)
hold on
errorbar(EXPDATA.time{2}, EXPDATA.mean{2}*k_ins_nM, EXPDATA.SD{2}*k_ins_nM,'Color', color11mM,...
    'LineStyle', 'n', 'Marker',marker,'MarkerFaceColor',color11mM, 'MarkerSize',8)
hold on
plot(simTime,plotData_11mM.Imeasured_opt,'Color', color11mM,'Linewidth',2);
hold off

xlim([0 simTime(end)])
xticks([0 48 96 144 192 240 288 336])
title('Measured, liver-islet co-culture','FontSize',10)
xlabel('Time (h)','FontSize',20);
ylabel('Insulin concentration (nM)','FontSize',20)
set(gca,'TickDir','out','FontSize',20);
box off

%% Simulate for single-liver culture, hyperglycemia (11 mM)

simTime=[0:0.01:time_end];
nTime = length(simTime);

param(:,parIndex.i_islets)=0;

% Only glucose offset at day 1 in the experiments when only liver is present

param(:,parIndex.i_delta_G_1)=param(:,parIndex.i_delta_G_1_liver);
param(:,parIndex.i_delta_G_7)=0;
param(:,parIndex.i_delta_G_13)=0;

param(:,parIndex.i_delta_I_1)=0;
param(:,parIndex.i_delta_I_7)=0;
param(:,parIndex.i_delta_I_13)=0;

% Simulate for only liver

plotData_singleliver = simulateCycles(param,parIndex,modelname,time_intervals,COSTOPTIONS,N_exchanges,k_ins_nM);


plotData_singleliver.Gmeasured_max(1:nTime)=max(plotData_singleliver.Gmeasured,[],1);
plotData_singleliver.Gmeasured_min(1:nTime)=min(plotData_singleliver.Gmeasured,[],1);
plotData_singleliver.Imeasured_max(1:nTime)=max(plotData_singleliver.Imeasured,[],1);
plotData_singleliver.Imeasured_min(1:nTime)=min(plotData_singleliver.Imeasured,[],1);

% Simulate for the optimal parameter values

% Change values for the optimal parameter

param_opt_onlyliver=param_opt;

param_opt_onlyliver(parIndex.i_islets)=0;
param_opt_onlyliver(parIndex.i_delta_G_1)=param_opt_onlyliver(parIndex.i_delta_G_1_liver);
param_opt_onlyliver(parIndex.i_delta_G_7)=0;
param_opt_onlyliver(parIndex.i_delta_G_13)=0;

param_opt_onlyliver(parIndex.i_delta_I_1)=0;
param_opt_onlyliver(parIndex.i_delta_I_7)=0;
param_opt_onlyliver(parIndex.i_delta_I_13)=0;

plotData_singleliver_opt = simulateCycles(param_opt_onlyliver,parIndex,modelname,time_intervals,COSTOPTIONS,N_exchanges,k_ins_nM);

plotData_singleliver.Gmeasured_opt=plotData_singleliver_opt.Gmeasured;

figure()

% Measured glucose
hold on

area1=[simTime,fliplr(simTime)];
area2=[plotData_singleliver.Gmeasured_max,fliplr(plotData_singleliver.Gmeasured_min)];
f=fill(area1,area2,color_onlyliver,'LineStyle','none');
alpha(f,.2)

hold on
errorbar(EXPDATA.time{3}, EXPDATA.mean{3}, EXPDATA.SD{3},'Color', color_onlyliver,...
    'LineStyle', 'n', 'Marker',marker,'MarkerFaceColor',color_onlyliver, 'MarkerSize',8)
hold on
plot(simTime,plotData_singleliver.Gmeasured_opt,'Color', color_onlyliver,'Linewidth',2);
hold off

xlim([0 simTime(end)])
ylim([0 12])
xticks([0 48 96 144 192 240 288 336])
title('Measured, single-liver culture')
xlabel('Time (h)','FontSize',20);
ylabel('Glucose concentration (mM)','FontSize',20)
set(gca,'TickDir','out','FontSize',20);
box off
