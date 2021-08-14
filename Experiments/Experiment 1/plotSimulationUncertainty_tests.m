% Plot model simulations with uncertainty based on the estimated parameter
% values

clear
clc
close all
format long
format compact

% Declare global variables

global EXPDATA_G
global EXPDATA_I
global EXPDATA_G_hep
global modelName
global icOrig0
global pNames
global param
global OPTIONS
global COSTOPTIONS
global EXPDATA
global SIMTIME      % Time axis for simulation
global OPTIME       % Time axis for optimization
global N_VAR        % Number of variables to optimize
global VAR_INDEX    % Indexes of variables to be optimized
global CUTOFF
global dataPoints
global FID
global Optparam
global parIndex

%% Load design parameters for the 2-OC system

% % Parameters to scale back to humans (alt)
% 
% V_hep_scaled = 3.40E-06*1e5;
% V_islets_scaled = 9.38E-08*1e5;
% Total_V_m= (V_hep_scaled+V_islets_scaled);
% % V_m_hep_scaled = Total_V_m/2;
% % V_m_islets_scaled = Total_V_m/2;
% V_m_hep_scaled = Total_V_m/2;
% %V_m_islets_scaled = Total_V_m/2;
% V_m_islets_scaled = Total_V_m/2;
% V_m_blood = 5;


MPSDesignParameters



%% Load and plot experimental data

loadExperimentalData

modelName ='M1_blood_compartment';

optModel = SBmodel(strcat(modelName,'.txt')); 
SBPDmakeMEXmodel(optModel,modelName); 
[pNames, param] = SBparameters(optModel);
icOrig0 = SBinitialconditions(optModel);

SIMTIME=[0:0.01:384];

% Parameter indexes

%% Parameter indexes

parIndex.i_V_m_hep = ismember(pNames,'V_m_hep');
parIndex.i_V_hep = ismember(pNames,'V_hep');
parIndex.i_V_m_islets = ismember(pNames,'V_m_islets');
parIndex.i_V_m_blood = ismember(pNames,'V_m_blood');
parIndex.i_Q = ismember(pNames,'Q');
parIndex.i_EGP_hep = ismember(pNames,'EGP_hep');
parIndex.i_U_ii_hep = ismember(pNames,'U_ii_hep');
parIndex.i_S_i = ismember(pNames,'S_i');
parIndex.i_Sigma = ismember(pNames,'Sigma');
parIndex.i_Alpha = ismember(pNames,'Alpha');
parIndex.i_CL_i_hep = ismember(pNames,'CL_i_hep');
parIndex.i_V_sample_hep = ismember(pNames,'V_sample_hep');
parIndex.i_V_sample_islets = ismember(pNames,'V_sample_islets');

% Disease progession parameters
parIndex.i_Vm = ismember(pNames,'Vm');
parIndex.i_km = ismember(pNames,'km');
parIndex.i_d0 = ismember(pNames,'d0');
parIndex.i_r1 = ismember(pNames,'r1');
parIndex.i_r2 = ismember(pNames,'r2');
parIndex.i_tao_slow = ismember(pNames,'tao_slow');
parIndex.i_G_healthy = ismember(pNames,'G_healthy');
parIndex.i_hep = ismember(pNames,'hep');
parIndex.i_islets = ismember(pNames,'islets');
parIndex.i_kv = ismember(pNames,'kv');
parIndex.i_Vm_f_hep = ismember(pNames,'Vm_f_hep');
parIndex.i_Km_f_hep = ismember(pNames,'Km_f_hep');
parIndex.i_Vm_f_islets = ismember(pNames,'Vm_f_islets');
parIndex.i_Km_f_islets = ismember(pNames,'Km_f_islets');

% Media change parameters
parIndex.i_G0 = ismember(pNames,'G0');
parIndex.i_I0 = ismember(pNames,'I0');
parIndex.i_G_GTT = ismember(pNames,'G_GTT');
parIndex.i_I_GTT = ismember(pNames,'I_GTT');

parIndex.i_delta_G_1 = ismember(pNames,'delta_G_1');
parIndex.i_delta_G_7 = ismember(pNames,'delta_G_7');
parIndex.i_delta_G_13 = ismember(pNames,'delta_G_13');

parIndex.i_delta_I_1 = ismember(pNames,'delta_I_1');
parIndex.i_delta_I_7 = ismember(pNames,'delta_I_7');
parIndex.i_delta_I_13 = ismember(pNames,'delta_I_13');

% Media change parameters for each concentration

% 11 mM

parIndex.i_delta_G_1_11mM = ismember(pNames,'delta_G_1_11mM');
parIndex.i_delta_G_7_11mM = ismember(pNames,'delta_G_7_11mM');
parIndex.i_delta_G_13_11mM = ismember(pNames,'delta_G_13_11mM');

parIndex.i_delta_I_1_11mM = ismember(pNames,'delta_I_1_11mM');
parIndex.i_delta_I_7_11mM = ismember(pNames,'delta_I_7_11mM');
parIndex.i_delta_I_13_11mM = ismember(pNames,'delta_I_13_11mM');

% 5.5 mM

parIndex.i_delta_G_1_5p5mM = ismember(pNames,'delta_G_1_5p5mM');
parIndex.i_delta_G_7_5p5mM = ismember(pNames,'delta_G_7_5p5mM');
parIndex.i_delta_G_13_5p5mM = ismember(pNames,'delta_G_13_5p5mM');

parIndex.i_delta_I_1_5p5mM = ismember(pNames,'delta_I_1_5p5mM');
parIndex.i_delta_I_7_5p5mM = ismember(pNames,'delta_I_7_5p5mM');
parIndex.i_delta_I_13_5p5mM = ismember(pNames,'delta_I_13_5p5mM');

% Only liver

parIndex.i_delta_G_1_liver = ismember(pNames,'delta_G_1_liver');

% iC=[(11+Optparam(parIndex.i_delta_G_1_11mM))*param(parIndex.i_V_m_hep) ...
%     (11+Optparam(parIndex.i_delta_G_1_11mM))*param(parIndex.i_V_m_islets)...
%    (0+Optparam(parIndex.i_delta_I_1_11mM))*param(parIndex.i_V_m_hep) (0+Optparam(parIndex.i_delta_I_1_11mM))*param(parIndex.i_V_m_islets) 0 0 ...
%    5.5 0.0000000088];

% SBPDsimulate options
OPTIONS.outputFunction ='';
OPTIONS.silent = 0;
OPTIONS.tempend =       0.1;                     % EndTemp
OPTIONS.tempfactor =    0.1;                     % tempfactor
OPTIONS.maxtime =       120;                     % Max-time
OPTIONS.tolx =          1e-10;                   % TolX
OPTIONS.tolfun =        1e-10;                   % Tolfun
OPTIONS.MaxRestartPoints =  0;                   % Number of parallel valleys which are searched through


% New scaling parameters to translate to humans

% Parameters when scaling keeping the volume of the pacreas

% V_hep_scaled = 3.54E-07;
% V_m_hep_scaled = 2.25E-07;
% V_m_islets_scaled = 2.25E-07;
% Q_scaled  = 8.04E-06;

% Parameters when scaling keeping the volume of the liver

% V_islets_scaled = 9.38E-08;
% V_m_hep_scaled = 2.88E-06;
% V_m_islets_scaled = 2.88E-06;
% Q_scaled  = 1.44E-04;
% V_hep_scaled = 3.40E-06;

% Parameters scaled back to humans

% V_islets_scaled = 9.38E-08*1e5;
% V_m_hep_scaled = 5.3/2;
% V_m_islets_scaled = 5.3/2;
% Q_scaled  = 4.5*60;
% V_hep_scaled = 3.40E-06*1e5;
% 
% % Parameters to scale back to humans (alt)
% 
V_hep_scaled = 3.40E-06*1e5;
V_islets_scaled = 9.38E-08*1e5;

Total_V_m= (V_hep_scaled+V_islets_scaled);
V_m_hep_scaled = Total_V_m/4;
V_m_islets_scaled = Total_V_m/4;

% V_m_hep_scaled=V_hep_scaled/2;
% V_m_islets_scaled=V_islets_scaled/2;
V_m_blood = 5;
Total_V_m=V_m_hep_scaled+V_m_islets_scaled;

% 
% %Recirculation time in the human body is about 5 minutes=5/60 hours
% r_time=5/60;
% Q_scaled=r_time/Total_V_m;
% V_hep_scaled = 3.40E-06*1e5;
% 
% Parameters to scale back to humans (alt 2)
% V_hep_scaled = 3.40E-06*1e5;
% V_m_hep_scaled = 2.5;
% V_m_islets_scaled = 2.5;
% V_islets_scaled = 8.8E-09*1e5;
% Total_V_m= (V_m_hep_scaled+V_m_islets_scaled);

%Recirculation time in the human body is about 5 minutes=5/60 hours
r_time=5/60;
Q_scaled=(Total_V_m+V_m_blood)/r_time;

% Calculate with residence time
% res_time_scaled=0.028; %(h)
% Q_scaled=V_m_hep_scaled/res_time_scaled;

% 
file = ['parameterValues M1 23-Jun-2020 18:36:16.dat'];

%Recirculation time in the human body is about 5 minutes=5/60 hours
% r_time=5/60;
% Q_scaled=(V_m_blood)/r_time;
% Q_scaled=Q_scaled/10;
% % 
% % Initial parameters
% V_hep_scaled = 3.40E-06;
% V_m_hep_scaled = 300E-06;
% V_m_islets_scaled = 300E-06;
% V_islets_scaled = 8.8E-09;
% Q_scaled = 2.96E-04;


tmpHold = load(file);


if isempty(tmpHold)
    disp('Empty parameter dataset');
else
    X = getSpreadedParam(tmpHold(:,2:end));
    X_opt = tmpHold(end,2:end);
    
    X=X(1,1:end);
    
    % Change to new scaling
    V_m_blood=V_m_blood;
    X(parIndex.i_V_m_hep)=V_m_hep_scaled;
    X(parIndex.i_V_m_islets)=V_m_islets_scaled;
    X(parIndex.i_V_hep)=V_hep_scaled;
    X(parIndex.i_Q)=Q_scaled;
    X(parIndex.i_CL_i_hep)=1.25*X(parIndex.i_CL_i_hep);
    X(parIndex.i_delta_G_1_11mM)=0;
    X(parIndex.i_delta_I_1_11mM)=0;
    X(parIndex.i_S_i)=X(parIndex.i_S_i);
    X(parIndex.i_U_ii_hep)=X(parIndex.i_S_i);
    
    
    
    
    X=[X V_m_blood];    % Add the volume in the blood compartment
    
    X_opt=[X_opt V_m_blood];
    
    % Change also the optimal parameter values
    
%     X_opt(parIndex.i_V_hep)=V_hep_scaled;
%     X_opt(parIndex.i_V_m_hep)=V_m_hep_scaled;
%     X_opt(parIndex.i_V_m_islets)=V_m_islets;
%     X_opt(parIndex.i_Q)=Q_scaled;
%     
    plotFunctionUncertainty(X,modelName,parIndex,X_opt,V_islets_scaled)

end




