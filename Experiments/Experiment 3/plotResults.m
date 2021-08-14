% Plot results for Donor 3

clear
clc
close all
format long
format compact

% Declare global variables

global EXPDATA_G
global EXPDATA_I
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
global statesIndex
global statesNames

%% Load design parameters for the 2-OC system

MPSDesignParameters

%% Load and plot experimental data

loadExperimentalData

% Define number of media exchanges (time intervals of 48 h) to simulate

N_exchanges=7;

modelName ='MPSmodel';

optModel = IQMmodel(strcat(modelName,'.txt')); 
IQMmakeMEXmodel(optModel,modelName); 
[pNames, param] = IQMparameters(optModel);
icOrig0 = IQMinitialconditions(optModel);
[statesNames,ODEs,initialConditions] = IQMstates(optModel);

SIMTIME=[0:1:336];

% Parameter indexes

% Parameter indexes

loadParameterIndexes

% Load state indexes
loadStateIndexes

file = ['ParameterValuesExperiment3.dat'];

tmpHold = load(file);

if isempty(tmpHold)
    disp('Empty parameter dataset');
else
    X = getSpreadedParam(tmpHold(:,2:end));
    X_opt = tmpHold(end,2:end);
    plotFunctionUncertainty(X,modelName,parIndex,X_opt,N_exchanges,statesIndex)
end
