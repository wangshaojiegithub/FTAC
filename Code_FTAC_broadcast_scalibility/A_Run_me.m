clear;clc;

load('Output_KNNnet_DWS1weight_S1S2_NodeFailRate.mat');
load('Output_KNNnet_Stableweight_S1S2_NodeFailRate.mat');
load('Output_NWnet_Stableweight_S1S2_NodeFailRate.mat');
load('Output_NWnet_DWS1weight_S1S2_NodeFailRate.mat');
load('Output_BAnet_Stableweight_S1S2_NodeFailRate.mat');
load('Output_BAnet_DWS1weight_S1S2_NodeFailRate.mat');

for Num=50:50:1000
    for ni= 1:20
        if ni==1
            eval(['temp_vec=Output_KNNnet_DWS1weight_S1S2_NodeFailRate.S1S2_NodeFailRate_of_',num2str(Num),'_nodes.NFR',num2str(ni)]);
        else
            eval(['tempt_value=Output_KNNnet_DWS1weight_S1S2_NodeFailRate.S1S2_NodeFailRate_of_',num2str(Num),'_nodes.NFR',num2str(ni)]);
            temp_vec=[temp_vec,tempt_value];
        end
    end    
    temp_vec_average=mean(temp_vec,2);
    if Num==50
        Y_ours_KNN_NodeFailRate=temp_vec_average;
    else
        Y_ours_KNN_NodeFailRate=[Y_ours_KNN_NodeFailRate,temp_vec_average];
    end
    temp_vec=[];
end

for Num=50:50:1000
    for ni= 1:20
        if ni==1
            eval(['temp_vec=Output_KNNnet_Stableweight_S1S2_NodeFailRate.S1S2_NodeFailRate_of_',num2str(Num),'_nodes.NFR',num2str(ni)]);
        else
            eval(['tempt_value=Output_KNNnet_Stableweight_S1S2_NodeFailRate.S1S2_NodeFailRate_of_',num2str(Num),'_nodes.NFR',num2str(ni)]);
            temp_vec=[temp_vec,tempt_value];
        end
    end    
    temp_vec_average=mean(temp_vec,2);
    if Num==50
        Y_others_KNN_NodeFailRate=temp_vec_average;
    else
        Y_others_KNN_NodeFailRate=[Y_others_KNN_NodeFailRate,temp_vec_average];
    end
    temp_vec=[];
end

for Num=50:50:1000
    for ni= 1:20
        if ni==1
            eval(['temp_vec=Output_NWnet_DWS1weight_S1S2_NodeFailRate.S1S2_NodeFailRate_of_',num2str(Num),'_nodes.NFR',num2str(ni)]);
        else
            eval(['tempt_value=Output_NWnet_DWS1weight_S1S2_NodeFailRate.S1S2_NodeFailRate_of_',num2str(Num),'_nodes.NFR',num2str(ni)]);
            temp_vec=[temp_vec,tempt_value];
        end
    end    
    temp_vec_average=mean(temp_vec,2);
    if Num==50
        Y_ours_NW_NodeFailRate=temp_vec_average;
    else
        Y_ours_NW_NodeFailRate=[Y_ours_NW_NodeFailRate,temp_vec_average];
    end
    temp_vec=[];
end

for Num=50:50:1000
    for ni= 1:20
        if ni==1
            eval(['temp_vec=Output_NWnet_Stableweight_S1S2_NodeFailRate.S1S2_NodeFailRate_of_',num2str(Num),'_nodes.NFR',num2str(ni)]);
        else
            eval(['tempt_value=Output_NWnet_Stableweight_S1S2_NodeFailRate.S1S2_NodeFailRate_of_',num2str(Num),'_nodes.NFR',num2str(ni)]);
            temp_vec=[temp_vec,tempt_value];
        end
    end    
    temp_vec_average=mean(temp_vec,2);
    if Num==50
        Y_others_NW_NodeFailRate=temp_vec_average;
    else
        Y_others_NW_NodeFailRate=[Y_others_NW_NodeFailRate,temp_vec_average];
    end
    temp_vec=[];
end

for Num=50:50:1000
    for ni= 1:20
        if ni==1
            eval(['temp_vec=Output_BAnet_DWS1weight_S1S2_NodeFailRate.S1S2_NodeFailRate_of_',num2str(Num),'_nodes.NFR',num2str(ni)]);
        else
            eval(['tempt_value=Output_BAnet_DWS1weight_S1S2_NodeFailRate.S1S2_NodeFailRate_of_',num2str(Num),'_nodes.NFR',num2str(ni)]);
            temp_vec=[temp_vec,tempt_value];
        end
    end    
    temp_vec_average=mean(temp_vec,2);
    if Num==50
        Y_ours_BA_NodeFailRate=temp_vec_average;
    else
        Y_ours_BA_NodeFailRate=[Y_ours_BA_NodeFailRate,temp_vec_average];
    end
    temp_vec=[];
end

for Num=50:50:1000
    for ni= 1:20
        if ni==1
            eval(['temp_vec=Output_BAnet_Stableweight_S1S2_NodeFailRate.S1S2_NodeFailRate_of_',num2str(Num),'_nodes.NFR',num2str(ni)]);
        else
            eval(['tempt_value=Output_BAnet_Stableweight_S1S2_NodeFailRate.S1S2_NodeFailRate_of_',num2str(Num),'_nodes.NFR',num2str(ni)]);
            temp_vec=[temp_vec,tempt_value];
        end
    end    
    temp_vec_average=mean(temp_vec,2);
    if Num==50
        Y_others_BA_NodeFailRate=temp_vec_average;
    else
        Y_others_BA_NodeFailRate=[Y_others_BA_NodeFailRate,temp_vec_average];
    end
    temp_vec=[];
end

plot_t=50:50:1000; 
figure(1);
p_NodeFailRate=plot(plot_t,Y_ours_KNN_NodeFailRate,'o-','color','m','LineWidth',0.7); 
hold on;grid on;
p_NodeFailRate=plot(plot_t,Y_ours_NW_NodeFailRate,'^-','color','g','LineWidth',0.7);
hold on;grid on;
p_NodeFailRate=plot(plot_t,Y_ours_BA_NodeFailRate,'s-','color','b','LineWidth',0.7);
hold on;grid on;

p_NodeFailRate=plot(plot_t,Y_others_KNN_NodeFailRate,'o:','color','m','LineWidth',1.2); 
hold on;grid on;
p_NodeFailRate=plot(plot_t,Y_others_NW_NodeFailRate,'^:','color','g','LineWidth',1.2);
hold on;grid on;
p_NodeFailRate=plot(plot_t,Y_others_BA_NodeFailRate,'s:','color','b','LineWidth',1.2);
hold on;grid on;
%     xlim([0 t_terminal-1]);
xlabel('Number of nodes (n)');    
ylabel('Node failure rate (%)');
legend('Our ALGO for KNNnet','Our ALGO for NWnet','Our ALGO for BAnet','[33] ALGO for KNNnet','[33] ALGO for NWnet','[33] ALGO for BAnet');
  


load('Output_KNNnet_Stableweight_S1_MaxTimestep.mat');
load('Output_KNNnet_DWS1weight_S1_MaxTimestep.mat');
load('Output_NWnet_Stableweight_S1_MaxTimestep.mat');
load('Output_NWnet_DWS1weight_S1_MaxTimestep.mat');
load('Output_BAnet_Stableweight_S1_MaxTimestep.mat');
load('Output_BAnet_DWS1weight_S1_MaxTimestep.mat');

for Num=50:50:1000
    for ni= 1:20
        if ni==1
            eval(['temp_vec=Output_KNNnet_DWS1weight_S1_MaxTimestep.S1_MaxTimestep_of_',num2str(Num),'_nodes.MT',num2str(ni)]);
        else
            eval(['tempt_value=Output_KNNnet_DWS1weight_S1_MaxTimestep.S1_MaxTimestep_of_',num2str(Num),'_nodes.MT',num2str(ni)]);
            temp_vec=[temp_vec,tempt_value];
        end
    end    
    temp_vec_average=mean(temp_vec,2);
    if Num==50
        Y_ours_KNN_MaxTimestep=temp_vec_average;
    else
        Y_ours_KNN_MaxTimestep=[Y_ours_KNN_MaxTimestep,temp_vec_average];
    end
    temp_vec=[];
end

for Num=50:50:1000
    for ni= 1:20
        if ni==1
            eval(['temp_vec=Output_KNNnet_Stableweight_S1_MaxTimestep.S1_MaxTimestep_of_',num2str(Num),'_nodes.MT',num2str(ni)]);
        else
            eval(['tempt_value=Output_KNNnet_Stableweight_S1_MaxTimestep.S1_MaxTimestep_of_',num2str(Num),'_nodes.MT',num2str(ni)]);
            temp_vec=[temp_vec,tempt_value];
        end
    end    
    temp_vec_average=mean(temp_vec,2);
    if Num==50
        Y_others_KNN_MaxTimestep=temp_vec_average;
    else
        Y_others_KNN_MaxTimestep=[Y_others_KNN_MaxTimestep,temp_vec_average];
    end
    temp_vec=[];
end

for Num=50:50:1000
    for ni= 1:20
        if ni==1
            eval(['temp_vec=Output_NWnet_DWS1weight_S1_MaxTimestep.S1_MaxTimestep_of_',num2str(Num),'_nodes.MT',num2str(ni)]);
        else
            eval(['tempt_value=Output_NWnet_DWS1weight_S1_MaxTimestep.S1_MaxTimestep_of_',num2str(Num),'_nodes.MT',num2str(ni)]);
            temp_vec=[temp_vec,tempt_value];
        end
    end    
    temp_vec_average=mean(temp_vec,2);
    if Num==50
        Y_ours_NW_MaxTimestep=temp_vec_average;
    else
        Y_ours_NW_MaxTimestep=[Y_ours_NW_MaxTimestep,temp_vec_average];
    end
    temp_vec=[];
end

for Num=50:50:1000
    for ni= 1:20
        if ni==1
            eval(['temp_vec=Output_NWnet_Stableweight_S1_MaxTimestep.S1_MaxTimestep_of_',num2str(Num),'_nodes.MT',num2str(ni)]);
        else
            eval(['tempt_value=Output_NWnet_Stableweight_S1_MaxTimestep.S1_MaxTimestep_of_',num2str(Num),'_nodes.MT',num2str(ni)]);
            temp_vec=[temp_vec,tempt_value];
        end
    end    
    temp_vec_average=mean(temp_vec,2);
    if Num==50
        Y_others_NW_MaxTimestep=temp_vec_average;
    else
        Y_others_NW_MaxTimestep=[Y_others_NW_MaxTimestep,temp_vec_average];
    end
    temp_vec=[];
end

for Num=50:50:1000
    for ni= 1:20
        if ni==1
            eval(['temp_vec=Output_BAnet_DWS1weight_S1_MaxTimestep.S1_MaxTimestep_of_',num2str(Num),'_nodes.MT',num2str(ni)]);
        else
            eval(['tempt_value=Output_BAnet_DWS1weight_S1_MaxTimestep.S1_MaxTimestep_of_',num2str(Num),'_nodes.MT',num2str(ni)]);
            temp_vec=[temp_vec,tempt_value];
        end
    end    
    temp_vec_average=mean(temp_vec,2);
    if Num==50
        Y_ours_BA_MaxTimestep=temp_vec_average;
    else
        Y_ours_BA_MaxTimestep=[Y_ours_BA_MaxTimestep,temp_vec_average];
    end
    temp_vec=[];
end

for Num=50:50:1000
    for ni= 1:20
        if ni==1
            eval(['temp_vec=Output_BAnet_Stableweight_S1_MaxTimestep.S1_MaxTimestep_of_',num2str(Num),'_nodes.MT',num2str(ni)]);
        else
            eval(['tempt_value=Output_BAnet_Stableweight_S1_MaxTimestep.S1_MaxTimestep_of_',num2str(Num),'_nodes.MT',num2str(ni)]);
            temp_vec=[temp_vec,tempt_value];
        end
    end    
    temp_vec_average=mean(temp_vec,2);
    if Num==50
        Y_others_BA_MaxTimestep=temp_vec_average;
    else
        Y_others_BA_MaxTimestep=[Y_others_BA_MaxTimestep,temp_vec_average];
    end
    temp_vec=[];
end

figure(2);
p_MaxTimestep=plot(plot_t,Y_ours_KNN_MaxTimestep,'o-','color','m','LineWidth',0.7); 
hold on;grid on;
p_MaxTimestep=plot(plot_t,Y_ours_NW_MaxTimestep,'^-','color','g','LineWidth',0.7);
hold on;grid on;
p_MaxTimestep=plot(plot_t,Y_ours_BA_MaxTimestep,'s-','color','b','LineWidth',0.7);
hold on;grid on;

p_MaxTimestep=plot(plot_t,Y_others_KNN_MaxTimestep,'o:','color','m','LineWidth',1.2); 
hold on;grid on;
p_MaxTimestep=plot(plot_t,Y_others_NW_MaxTimestep,'^:','color','g','LineWidth',1.2);
hold on;grid on;
p_MaxTimestep=plot(plot_t,Y_others_BA_MaxTimestep,'s:','color','b','LineWidth',1.2);
hold on;grid on;
%     xlim([0 t_terminal-1]);
xlabel('Number of nodes (n)');    
ylabel('Maximum time step (sec)');
legend('Our ALGO for KNNnet','Our ALGO for NWnet','Our ALGO for BAnet','[28] ALGO for KNNnet','[28] ALGO for NWnet','[28] ALGO for BAnet');
  

load('Output_KNNnet_Stableweight_S2_asyn_TerminalTime.mat');
load('Output_KNNnet_DWS1weight_S1S2_syn_TerminalTime.mat');
load('Output_NWnet_Stableweight_S2_asyn_TerminalTime.mat');
load('Output_NWnet_DWS1weight_S1S2_syn_TerminalTime.mat');
load('Output_BAnet_Stableweight_S2_asyn_TerminalTime.mat');
load('Output_BAnet_DWS1weight_S1S2_syn_TerminalTime.mat');

for Num=50:50:1000
    for ni= 1:20
        if ni==1
            eval(['temp_vec=Output_KNNnet_DWS1weight_S1S2_syn_TerminalTime.S1S2_syn_TerminalTime_of_',num2str(Num),'_nodes.TT',num2str(ni)]);
        else
            eval(['tempt_value=Output_KNNnet_DWS1weight_S1S2_syn_TerminalTime.S1S2_syn_TerminalTime_of_',num2str(Num),'_nodes.TT',num2str(ni)]);
            temp_vec=[temp_vec,tempt_value];
        end
    end    
    temp_vec_average=mean(temp_vec,2);
    if Num==50
        Y_ours_KNN_TerminalTime=temp_vec_average;
    else
        Y_ours_KNN_TerminalTime=[Y_ours_KNN_TerminalTime,temp_vec_average];
    end
    temp_vec=[];
end
Y_ours_KNN_scalability=Y_ours_KNN_TerminalTime/Y_ours_KNN_TerminalTime(1);

for Num=50:50:1000
    for ni= 1:20
        if ni==1
            eval(['temp_vec=Output_KNNnet_Stableweight_S2_asyn_TerminalTime.S2_asyn_TerminalTime_of_',num2str(Num),'_nodes.TT',num2str(ni)]);
        else
            eval(['tempt_value=Output_KNNnet_Stableweight_S2_asyn_TerminalTime.S2_asyn_TerminalTime_of_',num2str(Num),'_nodes.TT',num2str(ni)]);
            temp_vec=[temp_vec,tempt_value];
        end
    end    
    temp_vec_average=mean(temp_vec,2);
    if Num==50
        Y_others_KNN_TerminalTime=temp_vec_average;
    else
        Y_others_KNN_TerminalTime=[Y_others_KNN_TerminalTime,temp_vec_average];
    end
    temp_vec=[];
end
Y_others_KNN_scalability=Y_others_KNN_TerminalTime/Y_others_KNN_TerminalTime(1);

for Num=50:50:1000
    for ni= 1:20
        if ni==1
            eval(['temp_vec=Output_NWnet_DWS1weight_S1S2_syn_TerminalTime.S1S2_syn_TerminalTime_of_',num2str(Num),'_nodes.TT',num2str(ni)]);
        else
            eval(['tempt_value=Output_NWnet_DWS1weight_S1S2_syn_TerminalTime.S1S2_syn_TerminalTime_of_',num2str(Num),'_nodes.TT',num2str(ni)]);
            temp_vec=[temp_vec,tempt_value];
        end
    end    
    temp_vec_average=mean(temp_vec,2);
    if Num==50
        Y_ours_NW_TerminalTime=temp_vec_average;
    else
        Y_ours_NW_TerminalTime=[Y_ours_NW_TerminalTime,temp_vec_average];
    end
    temp_vec=[];
end
Y_ours_NW_scalability=Y_ours_NW_TerminalTime/Y_ours_NW_TerminalTime(1);

for Num=50:50:1000
    for ni= 1:20
        if ni==1
            eval(['temp_vec=Output_NWnet_Stableweight_S2_asyn_TerminalTime.S2_asyn_TerminalTime_of_',num2str(Num),'_nodes.TT',num2str(ni)]);
        else
            eval(['tempt_value=Output_NWnet_Stableweight_S2_asyn_TerminalTime.S2_asyn_TerminalTime_of_',num2str(Num),'_nodes.TT',num2str(ni)]);
            temp_vec=[temp_vec,tempt_value];
        end
    end    
    temp_vec_average=mean(temp_vec,2);
    if Num==50
        Y_others_NW_TerminalTime=temp_vec_average;
    else
        Y_others_NW_TerminalTime=[Y_others_NW_TerminalTime,temp_vec_average];
    end
    temp_vec=[];
end
Y_others_NW_scalability=Y_others_NW_TerminalTime/Y_others_NW_TerminalTime(1);

for Num=50:50:1000
    for ni= 1:20
        if ni==1
            eval(['temp_vec=Output_BAnet_DWS1weight_S1S2_syn_TerminalTime.S1S2_syn_TerminalTime_of_',num2str(Num),'_nodes.TT',num2str(ni)]);
        else
            eval(['tempt_value=Output_BAnet_DWS1weight_S1S2_syn_TerminalTime.S1S2_syn_TerminalTime_of_',num2str(Num),'_nodes.TT',num2str(ni)]);
            temp_vec=[temp_vec,tempt_value];
        end
    end    
    temp_vec_average=mean(temp_vec,2);
    if Num==50
        Y_ours_BA_TerminalTime=temp_vec_average;
    else
        Y_ours_BA_TerminalTime=[Y_ours_BA_TerminalTime,temp_vec_average];
    end
    temp_vec=[];
end
Y_ours_BA_scalability=Y_ours_BA_TerminalTime/Y_ours_BA_TerminalTime(1);

for Num=50:50:1000
    for ni= 1:20
        if ni==1
            eval(['temp_vec=Output_BAnet_Stableweight_S2_asyn_TerminalTime.S2_asyn_TerminalTime_of_',num2str(Num),'_nodes.TT',num2str(ni)]);
        else
            eval(['tempt_value=Output_BAnet_Stableweight_S2_asyn_TerminalTime.S2_asyn_TerminalTime_of_',num2str(Num),'_nodes.TT',num2str(ni)]);
            temp_vec=[temp_vec,tempt_value];
        end
    end    
    temp_vec_average=mean(temp_vec,2);
    if Num==50
        Y_others_BA_TerminalTime=temp_vec_average;
    else
        Y_others_BA_TerminalTime=[Y_others_BA_TerminalTime,temp_vec_average];
    end
    temp_vec=[];
end
Y_others_BA_scalability=Y_others_BA_TerminalTime/Y_others_BA_TerminalTime(1);

figure(3);
p_scalability=plot(plot_t,Y_ours_KNN_scalability,'o-','color','m','LineWidth',0.7); 
hold on;grid on;
p_scalability=plot(plot_t,Y_ours_NW_scalability,'^-','color','g','LineWidth',0.7);
hold on;grid on;
p_scalability=plot(plot_t,Y_ours_BA_scalability,'s-','color','b','LineWidth',0.7);
hold on;grid on;

p_scalability=plot(plot_t,Y_others_KNN_scalability,'o:','color','m','LineWidth',1.2); 
hold on;grid on;
p_scalability=plot(plot_t,Y_others_NW_scalability,'^:','color','g','LineWidth',1.2);
hold on;grid on;
p_scalability=plot(plot_t,Y_others_BA_scalability,'s:','color','b','LineWidth',1.2);
hold on;grid on;
    xlim([0 1000]);
xlabel('Number of nodes (n)');    
ylabel('Scalability');
legend('Our ALGO for KNNnet','Our ALGO for NWnet','Our ALGO for BAnet','[33] ALGO for KNNnet','[33] ALGO for NWnet','[33] ALGO for BAnet');
  

load('NWnet_AveLengthofLMLs.mat');
load('KNNnet_AveLengthofLMLs.mat');
load('BAnet_AveLengthofLMLs.mat');

for Num=50:50:1000
    for ni= 1:20
        if ni==1
            eval(['temp_vec=KNNnet_AveLengthofLMLs.KNNnet_AveLengthofLMLs_',num2str(Num),'_nodes.A',num2str(ni)]);
        else
            eval(['tempt_value=KNNnet_AveLengthofLMLs.KNNnet_AveLengthofLMLs_',num2str(Num),'_nodes.A',num2str(ni)]);
            temp_vec=[temp_vec,tempt_value];
        end
    end    
    temp_vec_average=mean(temp_vec,2);
    if Num==50
        Y_KNN_AveLengthofLCLs=temp_vec_average;
    else
        Y_KNN_AveLengthofLCLs=[Y_KNN_AveLengthofLCLs,temp_vec_average];
    end
    temp_vec=[];
end

for Num=50:50:1000
    for ni= 1:20
        if ni==1
            eval(['temp_vec=NWnet_AveLengthofLMLs.NWnet_AveLengthofLMLs_',num2str(Num),'_nodes.A',num2str(ni)]);
        else
            eval(['tempt_value=NWnet_AveLengthofLMLs.NWnet_AveLengthofLMLs_',num2str(Num),'_nodes.A',num2str(ni)]);
            temp_vec=[temp_vec,tempt_value];
        end
    end    
    temp_vec_average=mean(temp_vec,2);
    if Num==50
        Y_NW_AveLengthofLCLs=temp_vec_average;
    else
        Y_NW_AveLengthofLCLs=[Y_NW_AveLengthofLCLs,temp_vec_average];
    end
    temp_vec=[];
end

for Num=50:50:1000
    for ni= 1:20
        if ni==1
            eval(['temp_vec=BAnet_AveLengthofLMLs.BAnet_AveLengthofLMLs_',num2str(Num),'_nodes.A',num2str(ni)]);
        else
            eval(['tempt_value=BAnet_AveLengthofLMLs.BAnet_AveLengthofLMLs_',num2str(Num),'_nodes.A',num2str(ni)]);
            temp_vec=[temp_vec,tempt_value];
        end
    end    
    temp_vec_average=mean(temp_vec,2);
    if Num==50
        Y_BA_AveLengthofLCLs=temp_vec_average;
    else
        Y_BA_AveLengthofLCLs=[Y_BA_AveLengthofLCLs,temp_vec_average];
    end
    temp_vec=[];
end

figure(4);
p_AveLengthofLCLs=plot(plot_t,Y_KNN_AveLengthofLCLs,'o-','color','m','LineWidth',0.7); 
hold on;grid on;
p_AveLengthofLCLs=plot(plot_t,Y_NW_AveLengthofLCLs,'^-','color','g','LineWidth',0.7);
hold on;grid on;
p_AveLengthofLCLs=plot(plot_t,Y_BA_AveLengthofLCLs,'s-','color','b','LineWidth',0.7);
hold on;grid on;

%     xlim([0 1000]);
xlabel('Number of nodes (n)');    
ylabel('Average length of LCLs');
legend('KNNnet','NWnet','BAnet');
  

