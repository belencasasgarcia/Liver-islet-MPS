function plotFunctionUncertainty(param, model, parIndex, param_opt,N_exchanges,statesIndex)

time_end=335.99;
N_samples_G0=4; %Number of initial glucose concentrations for GTT at day 13
% under hypoglycemia. For plotting purposes is set to 4 (4 values between
% maximal and minimal concentrations in the experiment)

loadExperimentalData

COSTOPTIONS = [];
COSTOPTIONS.maxnumsteps = 1e8;
COSTOPTIONS.abstol = 1e-10;
COSTOPTIONS.reltol = 1e-10;
COSTOPTIONS.minstep = 0;
COSTOPTIONS.maxstep = inf;

%% Plotting settings

color11mM = [0.6350, 0.0780, 0.1840];           % Red
color5p5mM = [0.4660, 0.6740, 0.1880];          % Blue
color2p8mM = [0.9290, 0.6940, 0.1250];	        % Orange

marker='o';    

k_ins_nM=1/144; %Correction factor to change insulin units to nM

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

iC_11mM=[11*param(parIndex.i_V_m_liver) 11*param(parIndex.i_V_m_pancreas)...
   0*param(parIndex.i_V_m_liver) 0*param(parIndex.i_V_m_pancreas) 0 0 ...
   5.5 0.0000000088];%% Calculate maximal and minimal values of the predictions

%% Simulate for liver-islet co-culture, hyperglycemia (11 mM)

param(:,parIndex.i_delta_G_13)=0;

plotData_11mM = simulateCycles(param,parIndex,modelname,time_intervals,COSTOPTIONS,N_exchanges,k_ins_nM);

% Calculate maximal and minimal values of the predictions

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


%% Simulate for the optimal parameter values

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

%% Simulate for liver-islet co-culture, normoglycemia (5.5 mM)

param(:,parIndex.i_G0)=5.5;

plotData_5p5mM = simulateCycles(param,parIndex,modelname,time_intervals,COSTOPTIONS,N_exchanges,k_ins_nM);

% Calculate maximal and minimal values of the predictions

plotData_5p5mM.Gmeasured_max(1:nTime)=max(plotData_5p5mM.Gmeasured,[],1);
plotData_5p5mM.Gmeasured_min(1:nTime)=min(plotData_5p5mM.Gmeasured,[],1);
plotData_5p5mM.Imeasured_max(1:nTime)=max(plotData_5p5mM.Imeasured,[],1);
plotData_5p5mM.Imeasured_min(1:nTime)=min(plotData_5p5mM.Imeasured,[],1);
plotData_5p5mM.Gliver_max(1:nTime)=max(plotData_5p5mM.Gliver,[],1);
plotData_5p5mM.Gliver_min(1:nTime)=min(plotData_5p5mM.Gliver,[],1);
plotData_5p5mM.Gislets_max(1:nTime)=max(plotData_5p5mM.Gislets,[],1);
plotData_5p5mM.Gislets_min(1:nTime)=min(plotData_5p5mM.Gislets,[],1);
plotData_5p5mM.Iliver_max(1:nTime)=max(plotData_5p5mM.Iliver,[],1);
plotData_5p5mM.Iliver_min(1:nTime)=min(plotData_5p5mM.Iliver,[],1);
plotData_5p5mM.Iislets_max(1:nTime)=max(plotData_5p5mM.Iislets,[],1);
plotData_5p5mM.Iislets_min(1:nTime)=min(plotData_5p5mM.Iislets,[],1);
plotData_5p5mM.U_ii_min(1:nTime)=min(plotData_5p5mM.U_ii,[],1);
plotData_5p5mM.U_ii_max(1:nTime)=max(plotData_5p5mM.U_ii,[],1);
plotData_5p5mM.U_id_min(1:nTime)=min(plotData_5p5mM.U_id,[],1);
plotData_5p5mM.U_id_max(1:nTime)=max(plotData_5p5mM.U_id,[],1);
plotData_5p5mM.Glucose_int_max(1:nTime) = max(plotData_5p5mM.Glucose_int,[],1);
plotData_5p5mM.Glucose_int_min(1:nTime) = min(plotData_5p5mM.Glucose_int,[],1);
plotData_5p5mM.SI_max(1:nTime) = max(plotData_5p5mM.SI,[],1);
plotData_5p5mM.SI_min(1:nTime) = min(plotData_5p5mM.SI,[],1);
plotData_5p5mM.Gslow_max(1:nTime) = max(plotData_5p5mM.Gslow,[],1);
plotData_5p5mM.Gslow_min(1:nTime) = min(plotData_5p5mM.Gslow,[],1);
plotData_5p5mM.Vislets_max(1:nTime) = max(plotData_5p5mM.Vislets,[],1);
plotData_5p5mM.Vislets_min(1:nTime) = min(plotData_5p5mM.Vislets,[],1);
plotData_5p5mM.Secretion_max(1:nTime) = max(plotData_5p5mM.Secretion,[],1);
plotData_5p5mM.Secretion_min(1:nTime) = min(plotData_5p5mM.Secretion,[],1);


%% Simulate for the optimal parameter values

param_opt(:,parIndex.i_G0)=5.5;

plotData_5p5mM_opt = simulateCycles(param_opt,parIndex,modelname,time_intervals,COSTOPTIONS,N_exchanges,k_ins_nM);

plotData_5p5mM.Gmeasured_opt = plotData_5p5mM_opt.Gmeasured;
plotData_5p5mM.Imeasured_opt = plotData_5p5mM_opt.Imeasured;
plotData_5p5mM.Gliver_opt = plotData_5p5mM_opt.Gliver;
plotData_5p5mM.Gislets_opt = plotData_5p5mM_opt.Gislets;
plotData_5p5mM.Iliver_opt = plotData_5p5mM_opt.Iliver;
plotData_5p5mM.Iislets_opt = plotData_5p5mM_opt.Iislets;
plotData_5p5mM.U_ii_opt = plotData_5p5mM_opt.U_ii;
plotData_5p5mM.U_id_opt = plotData_5p5mM_opt.U_id;

plotData_5p5mM.Glucose_int_opt = plotData_5p5mM_opt.Glucose_int;
plotData_5p5mM.SI_opt = plotData_5p5mM_opt.SI;
plotData_5p5mM.Gslow_opt = plotData_5p5mM_opt.Gslow;
plotData_5p5mM.Vislets_opt = plotData_5p5mM_opt.Vislets;
plotData_5p5mM.Secretion_opt = plotData_5p5mM_opt.Secretion;

%% 2.8 mM

param(:,parIndex.i_G0)=2.8;
param_IC=repmat(param,N_samples_G0,1); %Sample N_samples values of glucose concentrations
% between (mean-SEM,mean+SEM)

% Simulate for concentration 2.8 mM

iC_2p8mM=[2.8*param(parIndex.i_V_m_liver) 2.8*param(parIndex.i_V_m_pancreas)...
   0*param(parIndex.i_V_m_liver) 0*param(parIndex.i_V_m_pancreas) 0 0 ...
   5.5 0.0000000088];%% Calculate maximal and minimal values of the predictions

param_IC0=[linspace(-0.8448,0,N_samples_G0/2) linspace(0,0.8448,N_samples_G0/2)];
param_IC0=repmat(param_IC0,size(param,1),1);
param_IC0=reshape(param_IC0,[],1);

param_IC(:,parIndex.i_delta_G_13)=param_IC0;
param_IC(:,parIndex.i_G_GTT)=repmat(11.112,size(param,1)*N_samples_G0,1);

plotData_2p8mM = simulateCycles(param_IC,parIndex,modelname,time_intervals,COSTOPTIONS,N_exchanges,k_ins_nM);

% Calculate maximal and minimal values of the predictions

plotData_2p8mM.Gmeasured_max(1:nTime)=max(plotData_2p8mM.Gmeasured,[],1);
plotData_2p8mM.Gmeasured_min(1:nTime)=min(plotData_2p8mM.Gmeasured,[],1);
plotData_2p8mM.Imeasured_max(1:nTime)=max(plotData_2p8mM.Imeasured,[],1);
plotData_2p8mM.Imeasured_min(1:nTime)=min(plotData_2p8mM.Imeasured,[],1);
plotData_2p8mM.Gliver_max(1:nTime)=max(plotData_2p8mM.Gliver,[],1);
plotData_2p8mM.Gliver_min(1:nTime)=min(plotData_2p8mM.Gliver,[],1);
plotData_2p8mM.Gislets_max(1:nTime)=max(plotData_2p8mM.Gislets,[],1);
plotData_2p8mM.Gislets_min(1:nTime)=min(plotData_2p8mM.Gislets,[],1);
plotData_2p8mM.Iliver_max(1:nTime)=max(plotData_2p8mM.Iliver,[],1);
plotData_2p8mM.Iliver_min(1:nTime)=min(plotData_2p8mM.Iliver,[],1);
plotData_2p8mM.Iislets_max(1:nTime)=max(plotData_2p8mM.Iislets,[],1);
plotData_2p8mM.Iislets_min(1:nTime)=min(plotData_2p8mM.Iislets,[],1);
plotData_2p8mM.U_ii_min(1:nTime)=min(plotData_2p8mM.U_ii,[],1);
plotData_2p8mM.U_ii_max(1:nTime)=max(plotData_2p8mM.U_ii,[],1);
plotData_2p8mM.U_id_min(1:nTime)=min(plotData_2p8mM.U_id,[],1);
plotData_2p8mM.U_id_max(1:nTime)=max(plotData_2p8mM.U_id,[],1);
plotData_2p8mM.Glucose_int_max(1:nTime) = max(plotData_2p8mM.Glucose_int,[],1);
plotData_2p8mM.Glucose_int_min(1:nTime) = min(plotData_2p8mM.Glucose_int,[],1);
plotData_2p8mM.SI_max(1:nTime) = max(plotData_2p8mM.SI,[],1);
plotData_2p8mM.SI_min(1:nTime) = min(plotData_2p8mM.SI,[],1);
plotData_2p8mM.Gslow_max(1:nTime) = max(plotData_2p8mM.Gslow,[],1);
plotData_2p8mM.Gslow_min(1:nTime) = min(plotData_2p8mM.Gslow,[],1);
plotData_2p8mM.Vislets_max(1:nTime) = max(plotData_2p8mM.Vislets,[],1);
plotData_2p8mM.Vislets_min(1:nTime) = min(plotData_2p8mM.Vislets,[],1);
plotData_2p8mM.Secretion_max(1:nTime) = max(plotData_2p8mM.Secretion,[],1);
plotData_2p8mM.Secretion_min(1:nTime) = min(plotData_2p8mM.Secretion,[],1);

% Save glucose and insulin predictions under hypoglycemia for posterior
% analyses.

%prediction_G_2p8mM=plotData_2p8mM.Gmeasured;
%prediction_I_2p8mM=plotData_2p8mM.Imeasured;

%save predictionHypoglycemia_G.mat prediction_G_2p8mM
%save predictionHypoglycemia_I.mat prediction_I_2p8mM

% Simulate for the optimal parameter value

param_opt(:,parIndex.i_G0)=2.8;

plotData_2p8mM_opt = simulateCycles(param_opt,parIndex,modelname,time_intervals,COSTOPTIONS,N_exchanges,k_ins_nM);

plotData_2p8mM.Gmeasured_opt = plotData_2p8mM_opt.Gmeasured;
plotData_2p8mM.Imeasured_opt = plotData_2p8mM_opt.Imeasured;
plotData_2p8mM.Gliver_opt = plotData_2p8mM_opt.Gliver;
plotData_2p8mM.Gislets_opt = plotData_2p8mM_opt.Gislets;
plotData_2p8mM.Iliver_opt = plotData_2p8mM_opt.Iliver;
plotData_2p8mM.Iislets_opt = plotData_2p8mM_opt.Iislets;
plotData_2p8mM.U_ii_opt = plotData_2p8mM_opt.U_ii;
plotData_2p8mM.U_id_opt = plotData_2p8mM_opt.U_id;

plotData_2p8mM.Glucose_int_opt = plotData_2p8mM_opt.Glucose_int;
plotData_2p8mM.SI_opt = plotData_2p8mM_opt.SI;
plotData_2p8mM.Gslow_opt = plotData_2p8mM_opt.Gslow;
plotData_2p8mM.Vislets_opt = plotData_2p8mM_opt.Vislets;
plotData_2p8mM.Secretion_opt = plotData_2p8mM_opt.Secretion;

% Glucose for 2.8 mM

figure()

area1=[simTime,fliplr(simTime)];
area2=[plotData_2p8mM.Gmeasured_max,fliplr(plotData_2p8mM.Gmeasured_min)];
f=fill(area1,area2,color2p8mM,'LineStyle','none');
alpha(f,.2)

hold on
errorbar(EXPDATA.time{7}, EXPDATA.mean{7}, EXPDATA.SD{7},'Color', color2p8mM,...
    'LineStyle', 'n', 'Marker',marker,'MarkerFaceColor',color2p8mM, 'MarkerSize',8)

hold on
plot(simTime,plotData_2p8mM.Gmeasured_opt,'Color', color2p8mM,'Linewidth',2);
hold off

xlim([287.9 335.9])
xlabel('Time (h)','FontSize',20);
ylabel('Glucose concentration (mM)','FontSize',20)
set(gca,'TickDir','out','FontSize',20);
box off

title('Glucose, hypoglycemia (prediction)','FontSize',14)

% Insulin for 2.8 mM

figure()
hold on

area1=[simTime,fliplr(simTime)];
area2=[plotData_2p8mM.Imeasured_max,fliplr(plotData_2p8mM.Imeasured_min)];
f=fill(area1,area2,color2p8mM,'LineStyle','none');
alpha(f,.2)

hold on

hold on
errorbar(EXPDATA.time{8}, EXPDATA.mean{8}*k_ins_nM, EXPDATA.SD{8}*k_ins_nM,'Color', color2p8mM,...
    'LineStyle', 'n', 'Marker',marker,'MarkerFaceColor',color2p8mM, 'MarkerSize',8)

hold on
plot(simTime,plotData_2p8mM.Imeasured_opt,'Color', color2p8mM,'Linewidth',2);
hold off

xlim([0 simTime(end)])
title('Measured')
xlabel('Time (h)','FontSize',20);
ylabel('Insulin concentration (nM)','FontSize',20)
set(gca,'TickDir','out','FontSize',20);
box off

title('Insulin, hypoglycemia (prediction)','FontSize',14)

%% Compare glucose consumption between the different concentrations

figure()

% 11 mM

area1=[simTime,fliplr(simTime)];
area2=[plotData_11mM.Gmeasured_max,fliplr(plotData_11mM.Gmeasured_min)];
f=fill(area1,area2,color11mM,'LineStyle','none');
alpha(f,.2)

hold on
errorbar(EXPDATA.time{1}, EXPDATA.mean{1}, EXPDATA.SD{1},'Color', color11mM,...
    'LineStyle', 'n', 'Marker',marker,'MarkerFaceColor',color11mM, 'MarkerSize',8)

hold on
plot(simTime,plotData_11mM.Gmeasured_opt,'Color', color11mM,'Linewidth',2);
hold on

% 5.5 mM

area1=[simTime,fliplr(simTime)];
area2=[plotData_5p5mM.Gmeasured_max,fliplr(plotData_5p5mM.Gmeasured_min)];
f=fill(area1,area2,color5p5mM,'LineStyle','none');
alpha(f,.2)

hold on
errorbar(EXPDATA.time{5}, EXPDATA.mean{5}, EXPDATA.SD{5},'Color', color5p5mM,...
    'LineStyle', 'n', 'Marker',marker,'MarkerFaceColor',color5p5mM, 'MarkerSize',8)
hold on

plot(simTime,plotData_5p5mM.Gmeasured_opt,'Color', color5p5mM,'Linewidth',2);

hold on

% % 2.8 mM 
% 
area1=[simTime,fliplr(simTime)];
area2=[plotData_2p8mM.Gmeasured_max,fliplr(plotData_2p8mM.Gmeasured_min)];
f=fill(area1,area2,color2p8mM,'LineStyle','none');
alpha(f,.2)

hold on
errorbar(EXPDATA.time{7}, EXPDATA.mean{7}, EXPDATA.SD{7},'Color', color2p8mM,...
    'LineStyle', 'n', 'Marker',marker,'MarkerFaceColor',color2p8mM, 'MarkerSize',8)

hold on
plot(simTime,plotData_2p8mM.Gmeasured_opt,'Color', color2p8mM,'Linewidth',2);
hold off

xlim([0 335.9])
xlabel('Time (h)','FontSize',20);
ylabel('Glucose concentration (mM)','FontSize',20)
set(gca,'TickDir','out','FontSize',20);
box off

title('Glucose, all glycemic regimes (calibration + prediction)','FontSize',14)

%% Compare insulin consumption between the different concentrations

figure()

% 11 mM

area1=[simTime,fliplr(simTime)];
area2=[plotData_11mM.Imeasured_max,fliplr(plotData_11mM.Imeasured_min)];
f=fill(area1,area2,color11mM,'LineStyle','none');
alpha(f,.2)

hold on
errorbar(EXPDATA.time{2}, EXPDATA.mean{2}*k_ins_nM, EXPDATA.SD{2}*k_ins_nM,'Color', color11mM,...
    'LineStyle', 'n', 'Marker',marker,'MarkerFaceColor',color11mM, 'MarkerSize',8)
hold on
plot(simTime,plotData_11mM.Imeasured_opt,'Color', color11mM,'Linewidth',2);

hold on

% 5.5 mM

area1=[simTime,fliplr(simTime)];
area2=[plotData_5p5mM.Imeasured_max,fliplr(plotData_5p5mM.Imeasured_min)];
f=fill(area1,area2,color5p5mM,'LineStyle','none');
alpha(f,.2)

hold on
errorbar(EXPDATA.time{6}, EXPDATA.mean{6}*k_ins_nM, EXPDATA.SD{6}*k_ins_nM,'Color', color5p5mM,...
    'LineStyle', 'n', 'Marker',marker,'MarkerFaceColor',color5p5mM, 'MarkerSize',8)
hold on
plot(simTime,plotData_5p5mM.Imeasured_opt,'Color', color5p5mM,'Linewidth',2);

hold on

% 2.8 mM 

area1=[simTime,fliplr(simTime)];
area2=[plotData_2p8mM.Imeasured_max,fliplr(plotData_2p8mM.Imeasured_min)];
f=fill(area1,area2,color2p8mM,'LineStyle','none');
alpha(f,.2)

hold on
errorbar(EXPDATA.time{8}, EXPDATA.mean{8}*k_ins_nM, EXPDATA.SD{8}*k_ins_nM,'Color', color2p8mM,...
    'LineStyle', 'n', 'Marker',marker,'MarkerFaceColor',color2p8mM, 'MarkerSize',8)

hold on
plot(simTime,plotData_2p8mM.Imeasured_opt,'Color', color2p8mM,'Linewidth',2);
hold off

xlim([0 335.9])
xlabel('Time (h)','FontSize',20);
ylabel('Insulin concentration (nM)','FontSize',20)
set(gca,'TickDir','out','FontSize',20);
box off

title('Glucose, all glycemic regimes (calibration + prediction)','FontSize',14)

%% Plot disease progression variables

%% 11mM
figure()

hold on
area1=[simTime,fliplr(simTime)];
area2=[plotData_11mM.SI_max,(fliplr(plotData_11mM.SI_min))];
f=fill(area1,area2,color11mM,'LineStyle','none','LineStyle','none');
alpha(f,.2)

hold on

plot(simTime,plotData_11mM.SI_opt,'Color',...
color11mM,'LineWidth',2);

hold on
area1=[simTime,fliplr(simTime)];
area2=[plotData_5p5mM.SI_max,(fliplr(plotData_5p5mM.SI_min))];
f=fill(area1,area2,color5p5mM,'LineStyle','none','LineStyle','none');
alpha(f,.2)

hold on

plot(simTime,plotData_5p5mM.SI_opt,'Color',...
color5p5mM,'LineWidth',2);

xlim([0 336])
title('Insulin sensitivity')
xlabel('Time (h)','FontSize',20);
ylabel('Insulin sensitivity (^{L}/_{mIU*h})','FontSize',20)
set(gca,'TickDir','out','FontSize',20);
box off

figure()

hold on
area1=[simTime,fliplr(simTime)];
area2=[plotData_11mM.Gslow_max,fliplr(plotData_11mM.Gslow_min)];
f=fill(area1,area2,color11mM,'LineStyle','none','LineStyle','none');
alpha(f,.2)

hold on

plot(simTime,plotData_11mM.Gslow_opt,'Color',...
color11mM,'LineWidth',2);

hold on
area1=[simTime,fliplr(simTime)];
area2=[plotData_5p5mM.Gslow_max,fliplr(plotData_5p5mM.Gslow_min)];
f=fill(area1,area2,color5p5mM,'LineStyle','none','LineStyle','none');
alpha(f,.2)

hold on

plot(simTime,plotData_5p5mM.Gslow_opt,'Color',...
color5p5mM,'LineWidth',2);

xlim([0 336])
title('Slow glucose')
xlabel('Time (h)','FontSize',20);
ylabel('Slow glucose (mM)','FontSize',20)
set(gca,'TickDir','out','FontSize',20);
box off

figure()

hold on
area1=[simTime,fliplr(simTime)];
area2=[plotData_11mM.Vislets_max*1e9,fliplr(plotData_11mM.Vislets_min*1e9)];
f=fill(area1,area2,color11mM,'LineStyle','none','LineStyle','none');
alpha(f,.2)

hold on

plot(simTime,plotData_11mM.Vislets_opt*1e9,'Color',...
color11mM,'LineWidth',2);

hold on
area1=[simTime,fliplr(simTime)];
area2=[plotData_5p5mM.Vislets_max*1e9,fliplr(plotData_5p5mM.Vislets_min*1e9)];
f=fill(area1,area2,color5p5mM,'LineStyle','none','LineStyle','none');
alpha(f,.2)

hold on

plot(simTime,plotData_5p5mM.Vislets_opt*1e9,'Color',...
color5p5mM,'LineWidth',2);

xlim([0 336])
title('Islets volume')
xlabel('Time (h)','FontSize',20);
ylabel('Islets volume (nL)','FontSize',20)
set(gca,'TickDir','out','FontSize',20);
box off

figure()

hold on
area1=[simTime,fliplr(simTime)];
area2=[100*plotData_11mM.Secretion_max,100*fliplr(plotData_11mM.Secretion_min)];
f=fill(area1,area2,color11mM,'LineStyle','none','LineStyle','none');
alpha(f,.2)

hold on
plot(simTime,100*plotData_11mM.Secretion_opt,'Color',...
color11mM,'LineWidth',2);

hold on
area1=[simTime,fliplr(simTime)];
area2=[100*plotData_5p5mM.Secretion_max,100*fliplr(plotData_5p5mM.Secretion_min)];
f=fill(area1,area2,color5p5mM,'LineStyle','none');
alpha(f,.2)

hold on
plot(simTime,100*plotData_5p5mM.Secretion_opt,'Color',...
color5p5mM,'LineWidth',2);

xlim([0 336])
title('Insulin secretion capacity')
xlabel('Time (h)','FontSize',20);
ylabel('Insulin secretion capacity (% per max)','FontSize',20);
set(gca,'TickDir','out','FontSize',20);
box off

