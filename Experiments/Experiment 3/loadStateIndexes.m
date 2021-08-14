% Define states indexes

statesIndex.i_NG_m_hep = ismember(statesNames,'NG_m_liver');  %Amount of glucose molecules in the liver compartment
statesIndex.i_NG_m_islets = ismember(statesNames,'NG_m_pancreas');  %Amount of glucose molecules in the pancreas compartment
statesIndex.i_NI_m_hep = ismember(statesNames,'NI_m_liver');  %Amount of glucose molecules in the liver compartment
statesIndex.i_NI_m_islets = ismember(statesNames,'NI_m_pancreas');  %Amount of glucose molecules in the pancreas compartment
statesIndex.time = ismember(statesNames,'time');  %Simulation time
statesIndex.Int_G = ismember(statesNames,'Int_G');  %Time average of glucose concentration
statesIndex.G_Slow = ismember(statesNames,'G_slow');  %Daily average of glucose concentration
statesIndex.V_islets = ismember(statesNames,'V_islets');  %Beta cell volume


