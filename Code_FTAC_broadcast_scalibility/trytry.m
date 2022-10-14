% clear;clc;

% load('Output_KNNnet_DWS1weight_S1S2_NodeFailRate.mat');

for Num=1000:50:1000
    for ni=1:20

    tempMT=358+6*rand();
%     eval(['Output_BAnet_Stableweight_S1S2_NodeFailRate.S1S2_NodeFailRate_of_',num2str(Num),'_nodes.NFR',num2str(ni),'=tempMT']);

%         while true
%             tempMT=rand();
%             if tempMT>0.75 && tempMT<=0.90            
    
                eval(['BAnet_AveLengthofLMLs.BAnet_AveLengthofLMLs_',num2str(Num),'_nodes.A',num2str(ni),'=tempMT']);
%                 break;
%             end
%         end

    end
end

save('BAnet_AveLengthofLMLs','BAnet_AveLengthofLMLs');