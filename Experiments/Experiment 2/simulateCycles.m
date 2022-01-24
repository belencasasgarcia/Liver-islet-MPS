function plotData = simulateCycles(param,parIndex,modelname,time_intervals,COSTOPTIONS,N_exchanges,k_ins_nM)
%Simulates model output for an experiment including N_exchanges media
%exchanges given the parameters param.

for j = 1 : size(param,1) % For each row of the param set.
    
    iC0=[(param(j,parIndex.i_G0)+param(j,parIndex.i_delta_G_1))*param(j,parIndex.i_V_m_liver) ...
            (param(j,parIndex.i_G0)+param(j,parIndex.i_delta_G_1))*param(j,parIndex.i_V_m_pancreas)...
            (param(j,parIndex.i_I0)+param(j,parIndex.i_delta_I_1))*param(j,parIndex.i_V_m_liver) ...
            (param(j,parIndex.i_I0)+param(j,parIndex.i_delta_I_1))*param(j,parIndex.i_V_m_pancreas) 0 0 ...
            5.5 0.0000000088 0 0];

    % Initial conditions for each media exchange to account for errors in
    % glucose and insulin measurements
        
    iC_G_I=repmat([param(j,parIndex.i_G0)*param(j,parIndex.i_V_m_liver) param(j,parIndex.i_G0)*param(j,parIndex.i_V_m_pancreas)...
        0 0],7,1);

    %Offset GTT day 7 
    iC_G_I(4,:)=iC_G_I(4,:)+[param(j,parIndex.i_delta_G_7)*param(j,parIndex.i_V_m_liver)...
    param(j,parIndex.i_delta_G_7)*param(j,parIndex.i_V_m_pancreas) ...
    param(j,parIndex.i_delta_I_7)*param(j,parIndex.i_V_m_liver)...
    param(j,parIndex.i_delta_I_7)*param(j,parIndex.i_V_m_pancreas)];

    %Offset GTT day 13
    iC_G_I(7,:)=iC_G_I(7,:)+[param(j,parIndex.i_delta_G_13)*param(j,parIndex.i_V_m_liver)...
    param(j,parIndex.i_delta_G_13)*param(j,parIndex.i_V_m_pancreas) ...
    param(j,parIndex.i_delta_I_13)*param(j,parIndex.i_V_m_liver)...
    param(j,parIndex.i_delta_I_13)*param(j,parIndex.i_V_m_pancreas)];
    
    try
        simData_11mM{1} = feval(char(modelname), [0:0.01:47.99], iC0, param(j,:), COSTOPTIONS);           
        %iC=zeros(N_exchanges,size(simData_11mM{1}.statevalues(end,:)));
        
        for i = 2 : N_exchanges
            iC(i,:)=simData_11mM{i-1}.statevalues(end,:); 
            iC(i,1:4)=iC_G_I(i,:);  %Include glucose and insulin offsets related to media exchanges
            simData_11mM{i} = feval(char(modelname), [0:0.01:47.99]+48*(i-1), iC(i,:), param(j,:), COSTOPTIONS);           
        end
        
    catch error
        disp(['Simulation crashed, @ simulation for data: ' num2str ' ... ->' ]);
        error = inf;
        disp(param);
        disp(error)
        return;
    end
    
    
    for i = 1 : N_exchanges
        plotData.Gmeasured(j,time_intervals(i,:))=simData_11mM{i}.variablevalues(:,1);
        plotData.Imeasured(j,time_intervals(i,:))=simData_11mM{i}.variablevalues(:,2)*k_ins_nM;
        plotData.Gliver(j,time_intervals(i,:))=simData_11mM{i}.variablevalues(:,3);
        plotData.Gislets(j,time_intervals(i,:))=simData_11mM{i}.variablevalues(:,4);
        plotData.Iliver(j,time_intervals(i,:))=simData_11mM{i}.variablevalues(:,5)*k_ins_nM;
        plotData.Iislets(j,time_intervals(i,:))=simData_11mM{i}.variablevalues(:,6)*k_ins_nM;
        plotData.U_ii(j,time_intervals(i,:))=simData_11mM{i}.variablevalues(:,27);
        plotData.U_id(j,time_intervals(i,:))=simData_11mM{i}.variablevalues(:,28);
        plotData.Glucose_int(j,time_intervals(i,:))=simData_11mM{i}.variablevalues(:,7);
        plotData.SI(j,time_intervals(i,:))=(simData_11mM{i}.variablevalues(:,8))/k_ins_nM;
        plotData.Gslow(j,time_intervals(i,:))=simData_11mM{i}.variablevalues(:,9);
        plotData.Vislets(j,time_intervals(i,:))=simData_11mM{i}.variablevalues(:,12);  
        plotData.Secretion(j,time_intervals(i,:))=simData_11mM{i}.variablevalues(:,14);  
        plotData.Uptake_ii(j,time_intervals(i,:))=simData_11mM{i}.variablevalues(:,29);  
        plotData.Uptake_id(j,time_intervals(i,:))=simData_11mM{i}.variablevalues(:,30);  
    end
    
    
end

end

