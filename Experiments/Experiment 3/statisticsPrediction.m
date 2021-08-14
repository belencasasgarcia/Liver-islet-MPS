% Calculate statistics for prediction, optimal prediction

clear
clc
close all

loadExperimentalData

load predictionHypoglycemia_G
load predictionHypoglycemia_I

EXPDATA=[];
L_GTT=48;
N_exchanges=7;
% Correct units
k_ins_nM=1/144; %Change insulin units to nM

marker='o';         


color2p8mM = [0.9290, 0.6940, 0.1250];	    % Yellow

simTime=[0:0.01:(L_GTT*N_exchanges-0.01)];

%% Glucose and insulin data 2.8 mM

% Glucose data 2.8 mM 
EXPDATA.time{1} = EXPDATA_G_2p8mM([1 3 4 5],1);    % Time (absolute) (hours)
EXPDATA.mean{1} = EXPDATA_G_2p8mM([1 3 4 5],2);    % Time (absolute) (hours)
EXPDATA.SD{1} = EXPDATA_G_2p8mM([1 3 4 5],4);      % Time (absolute) (hours)

% Insulin data 2.8 mM 
EXPDATA.time{2} = EXPDATA_I_2p8mM([3 4 5],1);    % Time (absolute) (hours)
EXPDATA.mean{2} = EXPDATA_I_2p8mM([3 4 5],2);    % Time (absolute) (hours)
EXPDATA.SD{2} = EXPDATA_I_2p8mM([3 4 5],4);      % Time (absolute) (hours)

% Calculate number of data points for statistical tests
dataPoints=0;

OPTVAR_INDEX=[1 2];

for i=1:size(OPTVAR_INDEX,2)
    dataPoints=dataPoints+length(EXPDATA.time{i});
end

CUTOFF = chi2inv(0.95,dataPoints);

% Calculate cost value

tmpError=zeros(1,size(prediction_G_2p8mM,1));
tmpError_G=zeros(1,size(prediction_G_2p8mM,1));
tmpError_I=zeros(1,size(prediction_I_2p8mM,1));

for i=1:1:size(tmpError,2)
    % Glucose prediction
    index_int1=ismembertol(simTime,EXPDATA.time{1});    
    tmpError_G(i) = sum(((prediction_G_2p8mM(i,index_int1)'- EXPDATA.mean{1}).^2)./(EXPDATA.SD{1}).^2);
    index_int2=ismembertol(simTime,EXPDATA.time{2}); 
    tmpError_I(i) = sum(((prediction_I_2p8mM(i,index_int2)'- EXPDATA.mean{2}).^2)./(EXPDATA.SD{2}).^2);
    tmpError(i)=tmpError_G(i)+tmpError_I(i);
end

figure() 

for j = 1 : size(tmpError,2)
    plot(simTime,prediction_G_2p8mM(j,:),'Color',color2p8mM)
    hold on
end

hold on

errorbar(EXPDATA.time{1}, EXPDATA.mean{1}, EXPDATA.SD{1},'Color', color2p8mM,...
    'LineStyle', 'n', 'Marker',marker,'MarkerFaceColor',color2p8mM, 'MarkerSize',8)
title('Insulin, liver-islet co-culture, hypoglycemia', 'FontSize',10)
xlabel('Time (h)');
ylabel('Glucose concentration (mM)','FontSize',20)
xticks([0 48 96 144 192 240 288 336])

figure()

for j = 1 : size(tmpError,2)
    plot(simTime,prediction_I_2p8mM(j,:),'Color',color2p8mM)
    hold on
end

hold on

errorbar(EXPDATA.time{2}, EXPDATA.mean{2}*k_ins_nM, EXPDATA.SD{2}*k_ins_nM,'Color', color2p8mM,...
    'LineStyle', 'n', 'Marker',marker,'MarkerFaceColor',color2p8mM, 'MarkerSize',8)
title('Insulin, liver-islet co-culture, hypoglycemia', 'FontSize',10)
xlabel('Time (h)');
ylabel('Insulin concentration (nM)','FontSize',20)
xticks([0 48 96 144 192 240 288 336])



figure()

for j = 1 : size(tmpError,2)
    plot(simTime,prediction_I_2p8mM(j,:),'Color',color2p8mM)
    hold on
end

hold on

errorbar(EXPDATA.time{2}, EXPDATA.mean{2}*k_ins_nM, EXPDATA.SD{2}*k_ins_nM,'Color', color2p8mM,...
    'LineStyle', 'n', 'Marker',marker,'MarkerFaceColor',color2p8mM, 'MarkerSize',8)
title('Insulin, liver-islet co-culture, hypoglycemia', 'FontSize',10)
xlabel('Time (h)');
ylabel('Insulin concentration (nM)','FontSize',20)
xticks([0 48 96 144 192 240 288 336])

[min_error, min_index]=min(tmpError);    %Index of the optimal prediction

