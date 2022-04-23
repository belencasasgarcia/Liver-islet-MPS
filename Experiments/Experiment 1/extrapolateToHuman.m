function [parHuman V_islets_scaled] = extrapolateToHuman(parMPS,V_hep,V_islets,parIndex,adjustSigma,adjustCL)

% Extrapolate parameter values in the MPS (parMPS) to human values
% (parHuman)

% Define parameter values

V_hep_scaled = V_hep*1e5;             % Upscale based on downscaling factor (100000)
V_islets_scaled = V_islets*1e5;       % Upscale based on downscaling factor (100000)
V_m_plasma=0.58*5.1;                  % Blood volume in humans is approximately 5.1 L with a plasma proportion of 58%
V_m_liver_scaled = V_m_plasma/2;      % Half of plasma volume in each culture compartment  
V_m_pancreas_scaled = V_m_plasma/2;
Total_V_m=V_m_liver_scaled+V_m_pancreas_scaled;

%Recirculation time in the human body is about 5 minutes=5/60 hours
r_time=5/60;
Q_scaled=(Total_V_m)/r_time;

parHuman=parMPS;
parHuman(:,parIndex.i_V_m_liver)=V_m_liver_scaled;
parHuman(:,parIndex.i_V_m_pancreas)=V_m_pancreas_scaled;
parHuman(:,parIndex.i_V_hep)=V_hep_scaled;
parHuman(:,parIndex.i_Q)=Q_scaled;
parHuman(:,parIndex.i_CL_i_hep)=parMPS(:,parIndex.i_CL_i_hep)*2;        % Assuming that the liver stands for approximately 50% of total insulin clearance
parHuman(:,parIndex.i_delta_G_1_11mM)=0;
parHuman(:,parIndex.i_delta_I_1_11mM)=0;
parHuman(:,parIndex.i_S_i)=parMPS(:,parIndex.i_S_i)*2.22;               % Assuming that the liver is responsible for approximately 45% of the total postprandial glucose uptake in humans
parHuman(:,parIndex.i_EG0_hep)=parMPS(:,parIndex.i_EG0_hep)*2.22;

if (adjustSigma)
    parHuman(:,parIndex.i_Sigma)=parHuman(:,parIndex.i_Sigma)/2;
end

if (adjustCL)
    parHuman(:,parIndex.i_CL_i_hep)=parHuman(:,parIndex.i_CL_i_hep)*4.23;

end

