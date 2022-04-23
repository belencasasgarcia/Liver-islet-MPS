function plotFunctionUncertainty(param, model, parIndex, param_opt,V_islets_scaled)

loadExperimentalData

plotExperimentalData(EXPDATA)

COSTOPTIONS = [];
COSTOPTIONS.maxnumsteps = 1e8;
COSTOPTIONS.abstol = 1e-10;
COSTOPTIONS.reltol = 1e-10;
COSTOPTIONS.minstep = 0;
COSTOPTIONS.maxstep = inf;

time_end=4;

k_conv_I=1/144;

G_start_liver=11;
G_start_islets=11;

%k_conv_I=1;

%% Plotting settings

color11mM = [0.6350, 0.0780, 0.1840];           % Red
color_onlyliver = [0, 0.4470, 0.7410];          % Blue
color_onlyislets = [0.8500, 0.3250, 0.0980];	% Orange
color_DallaMan = [0 0 0];                       % Black


marker='o';              % Same marker for all concentrations

% 11 mM glucose

modelname{1} = char(model);
simTime=[0:0.01:time_end];
nTime = length(simTime);

% Simulate for concentration 11 mM


for j = 1 : size(param,1) % For each row of the param set.
    
    try            
            %Initial conditions
            
            iC_11mM=[(G_start_liver+param(j,parIndex.i_delta_G_1_11mM))*param(j,parIndex.i_V_m_liver) ...
            (G_start_islets+param(j,parIndex.i_delta_G_1_11mM))*param(j,parIndex.i_V_m_pancreas)...
            (0+param(j,parIndex.i_delta_I_1_11mM))*param(j,parIndex.i_V_m_pancreas) ...
            (0+param(j,parIndex.i_delta_I_1_11mM))*param(j,parIndex.i_V_m_pancreas) 0 0 ...
            5.5 V_islets_scaled];
        
            simData_11mM = feval(char(modelname), simTime, iC_11mM, param(j,:), COSTOPTIONS);        
        
    catch error
        disp(['Simulation crashed, @ simulation for data: ' num2str(j) ' ... ->' ]);
        error = inf;
        disp(param);
        disp(error)
        return;
    end
    
    plotData.Gmeasured_11mM(j,1:nTime) = simData_11mM.variablevalues(:,1);
    plotData.Imeasured_11mM(j,1:nTime) = simData_11mM.variablevalues(:,2);
    plotData.Gliver_11mM(j,1:nTime) = simData_11mM.variablevalues(:,3);
    plotData.Gislets_11mM(j,1:nTime) = simData_11mM.variablevalues(:,4);
    plotData.Iliver_11mM(j,1:nTime) = simData_11mM.variablevalues(:,5);
    plotData.Iislets_11mM(j,1:nTime) = simData_11mM.variablevalues(:,6);
    plotData.U_ii_11mM(j,1:nTime) = simData_11mM.variablevalues(:,27);
    plotData.U_id_11mM(j,1:nTime) = simData_11mM.variablevalues(:,28);
    
end

%% Calculate maximal and minimal values of the predictions

plotData.Gmeasured_11mM_max(1:nTime)=max(plotData.Gmeasured_11mM,[],1);
plotData.Gmeasured_11mM_min(1:nTime)=min(plotData.Gmeasured_11mM,[],1);
plotData.Imeasured_11mM_max(1:nTime)=max(plotData.Imeasured_11mM,[],1);
plotData.Imeasured_11mM_min(1:nTime)=min(plotData.Imeasured_11mM,[],1);
plotData.Gliver_11mM_max(1:nTime)=max(plotData.Gliver_11mM,[],1);
plotData.Gliver_11mM_min(1:nTime)=min(plotData.Gliver_11mM,[],1);
plotData.Gislets_11mM_max(1:nTime)=max(plotData.Gislets_11mM,[],1);
plotData.Gislets_11mM_min(1:nTime)=min(plotData.Gislets_11mM,[],1);
plotData.Iliver_11mM_max(1:nTime)=max(plotData.Iliver_11mM,[],1);
plotData.Iliver_11mM_min(1:nTime)=min(plotData.Iliver_11mM,[],1);
plotData.Iislets_11mM_max(1:nTime)=max(plotData.Iislets_11mM,[],1);
plotData.Iislets_11mM_min(1:nTime)=min(plotData.Iislets_11mM,[],1);
plotData.U_ii_11mM_min(1:nTime)=min(plotData.U_ii_11mM,[],1);
plotData.U_ii_11mM_max(1:nTime)=max(plotData.U_ii_11mM,[],1);
plotData.U_id_11mM_min(1:nTime)=min(plotData.U_id_11mM,[],1);
plotData.U_id_11mM_max(1:nTime)=max(plotData.U_id_11mM,[],1);

%% Simulations fitted to experimental measurements

figure()

% Glucose in liver compartment
hold on
area1=[simTime,fliplr(simTime)];
area2=[plotData.Gliver_11mM_max,fliplr(plotData.Gliver_11mM_min)];
f=fill(area1,area2,color_onlyliver,'LineStyle','none');
alpha(f,.4)

hold off
%title('Liver compartment, after translation to human','FontSize',20)
xlabel('Time (h)');
ylabel('Glucose concentration (mM)','FontSize',20)
xlim([0 simTime(end)])
ylim([0 12])
set(gca,'TickDir','out','FontSize',20);
box off

% 11 mM

% Glucose in the islets compartment
figure()
hold on
area1=[simTime,fliplr(simTime)];
area2=[plotData.Gislets_11mM_max,fliplr(plotData.Gislets_11mM_min)];
f=fill(area1,area2,color_onlyislets,'LineStyle','none');
alpha(f,.4)

hold off
%title('Pancreas compartment, after translation to human','FontSize',20)
xlabel('Time (h)');
ylabel('Glucose concentration (mM)','FontSize',20)
xlim([0 simTime(end)])
ylim([0 12])
set(gca,'TickDir','out','FontSize',20);
box off

% Measured glucose
figure()
hold on

area1=[simTime,fliplr(simTime)];
area2=[plotData.Gmeasured_11mM_max,fliplr(plotData.Gmeasured_11mM_min)];
f=fill(area1,area2,color11mM,'LineStyle','none');
alpha(f,.4)

hold on
errorbar(G_mean_Dallaman(:,1), G_mean_Dallaman(:,2), SEM_G,'Color', color_DallaMan,...
    'LineStyle', 'n', 'Marker',marker,'MarkerFaceColor',color_DallaMan, 'MarkerSize',6)

xlim([0 simTime(end)])
ylim([0 12])
%title('Measured, after translation to human','FontSize',20)
xlabel('Time (h)','FontSize',20);
ylabel('Glucose concentration (mM)','FontSize',20)
set(gca,'TickDir','out','FontSize',20);
box off


% Insulin in the liver compartment
figure()
hold on
area1=[simTime,fliplr(simTime)];
area2=[plotData.Iliver_11mM_max*k_conv_I,fliplr(plotData.Iliver_11mM_min*k_conv_I)];
f=fill(area1,area2,color_onlyliver,'LineStyle','none');
alpha(f,.4)

hold off
xlabel('Time (h)');
ylabel('Insulin concentration (nM)','FontSize',20)
%title('Liver compartment, after translation to human','FontSize',20)
xlim([0 simTime(end)])
set(gca,'TickDir','out','FontSize',20);
box off

% Insulin in the islets compartment
figure()

hold on
area1=[simTime,fliplr(simTime)];
area2=[plotData.Iislets_11mM_max*k_conv_I,fliplr(plotData.Iislets_11mM_min*k_conv_I)];
f=fill(area1,area2,color_onlyislets,'LineStyle','none');
alpha(f,.4)

hold off
xlabel('Time (h)');
ylabel('Insulin concentration (nM)','FontSize',20)
%title('Pancreas compartment, after translation to human','FontSize',20)
xlim([0 simTime(end)])
set(gca,'TickDir','out','FontSize',20);
box off

% Measured insulin
figure()

hold on

area1=[simTime,fliplr(simTime)];
area2=[plotData.Imeasured_11mM_max*k_conv_I,fliplr(plotData.Imeasured_11mM_min*k_conv_I)];
f=fill(area1,area2,color11mM,'LineStyle','none');
alpha(f,.4)
hold on

errorbar(I_mean_Dallaman(:,1), I_mean_Dallaman(:,2)/1000, SEM_I/1000,'Color', color_DallaMan,...
    'LineStyle', 'n', 'Marker',marker,'MarkerFaceColor',color_DallaMan, 'MarkerSize',6)
% hold on
% plot(simTime,plotData.Imeasured_11mM_opt,'Color', color11mM,'Linewidth',2);

hold off

xlim([0 simTime(end)])
xlabel('Time (h)','FontSize',20);
ylabel('Insulin concentration (nM)','FontSize',20)
%title('Measured, after translation to human','FontSize',20)
set(gca,'TickDir','out','FontSize',20);
box off



