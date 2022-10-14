clear;clc;
% **************Stage 1：计算行和相等的 W=L 的左特征向量************
    
load('KNN_matrices.mat');    DatasetNumber=2;

for Num=50:50:500
    for ni=1:20
        disp('Num=');disp(Num);
        disp('ni=');disp(ni);
        
        if DatasetNumber==2
            eval(['A_temp=KNN_matrices.KNN_matrices_of_',num2str(Num),'_nodes.A',num2str(ni)]);
        end
        Din=diag(sum(A_temp'));     %  A_temp的对应的入度矩阵         
        W=Din-A_temp;               % 将A_temp转换成Laplacian矩阵(行和为0)
        
        
        %****** 将A_temp转换成文献[28]中的行随机矩阵(用于对比Stage 1的左特征向量计算结果)****
            A_weight=A_temp;
            for i=1:size(A_temp,1)
                A_weight(i,A_weight(i,:)~=0) = 1/Din(i,i);
            end
            W=0.5*eye(size(A_temp,1))+0.5*A_weight;
        %***************************************************************************
                
        
        %******************************************   接下来是主程序    ********************************************************

            miu=sum(W(1,:));                % 求行和
            I=eye(size(W,1));
            for vi=1:size(W,1)  
                eval(['eT',num2str(vi),'=I(vi,:)']);
            end 

        % ************************************* Stage 1 计算左特征向量 ************************************************************
            S1_NodeFailRate=0;  %  该变量用于记录Stage 1的失败节点率
            S1_MaxTimestep=0;   %  该变量用于记录Stage 1的最大时间步
            for  vi=1:size(W,1)    
                eval(['eTi=eT',num2str(vi)]);
                [di_real,ri_max,ri_max_source_vec,LeftVecUni_real_at_vi,LeftVecUni_at_vi,MaxTimestep]=computing_LeftEigenvector_from_channels(vi,eTi,W);

                % 以列向量形式表达:S1_di_vec和S1_ri_max_vec
                if vi==1
                    S1_di_vec=di_real;
                    S1_ri_max_vec=ri_max;
                    S1_MaxTimestep=MaxTimestep;   % 不同节点的最大时间步 取最大值
                else
                    S1_di_vec=[S1_di_vec;di_real];
                    S1_ri_max_vec=[S1_ri_max_vec;ri_max];
                    S1_MaxTimestep=max(S1_MaxTimestep,MaxTimestep);   % 不同节点的最大时间步 取最大值
                end    
                    % ****************** 测试在节点vi处,对于左特征向量x的计算是否正确 ***********
                        if miu==1
                            % 若W是行随机矩阵,失败包含两种情况:1.出现计算错误的节点  2. vi处的di与某个通道的rj不同(即初始条件p0属于文献中描述的零测集)
                            standard_result=roundn(LeftVecUni_real_at_vi,-2);
                            computed_result=roundn(LeftVecUni_at_vi,-2);
                            if ~isempty(find((standard_result-computed_result)~=0)) || length(find(ri_max_source_vec~=di_real))>0
                                S1_NodeFailRate=S1_NodeFailRate+1;
                            end
                            if vi==size(W,1)
                                S1_NodeFailRate=S1_NodeFailRate/size(W,1); % 计算节点失败率
                            end
                        else
                            % 若W不是行随机矩阵,仅统计计算出现错误的节点比例
                            standard_result=roundn(LeftVecUni_real_at_vi,-2);
                            computed_result=roundn(LeftVecUni_at_vi,-2);
                            if ~isempty(find((standard_result-computed_result)~=0))                
                                S1_NodeFailRate=S1_NodeFailRate+1;
                            end
                            if vi==size(W,1)
                                S1_NodeFailRate=S1_NodeFailRate/size(W,1); % 计算节点失败率
                            end
                        end
                    % *************************************************************************   
            end

            if DatasetNumber==2
                %         存储最大时间步数据,用于对比: 本文的Stage 1的算法和[28]的左特征向量算法,
                if miu==1
                    eval(['Output_KNNnet_Stableweight_S1_MaxTimestep.S1_MaxTimestep_of_',num2str(Num),'_nodes.MT',num2str(ni),'=S1_MaxTimestep']);
                    eval(['Output_KNNnet_Stableweight_S1_NodeFailRate.S1_NodeFailRate_of_',num2str(Num),'_nodes.NFR',num2str(ni),'=S1_NodeFailRate']);
                    save('Output_KNNnet_Stableweight_S1_MaxTimestep','Output_KNNnet_Stableweight_S1_MaxTimestep');
                    save('Output_KNNnet_Stableweight_S1_NodeFailRate','Output_KNNnet_Stableweight_S1_NodeFailRate');  
                else
                    eval(['Output_KNNnet_DWS1weight_S1_MaxTimestep.S1_MaxTimestep_of_',num2str(Num),'_nodes.MT',num2str(ni),'=S1_MaxTimestep']);
                    eval(['Output_KNNnet_DWS1weight_S1_NodeFailRate.S1_NodeFailRate_of_',num2str(Num),'_nodes.NFR',num2str(ni),'=S1_NodeFailRate']);
                    save('Output_KNNnet_DWS1weight_S1_MaxTimestep','Output_KNNnet_DWS1weight_S1_MaxTimestep');
                    save('Output_KNNnet_DWS1weight_S1_NodeFailRate','Output_KNNnet_DWS1weight_S1_NodeFailRate');  
                end
            end

            clearvars -except KNN_matrices DatasetNumber Num ni W LeftVecUni_real_at_vi S1_di_vec S1_ri_max_vec Output_KNNnet_DWS1weight_S1_MaxTimestep Output_KNNnet_DWS1weight_S1_NodeFailRate S1_MaxTimestep S1_NodeFailRate
          
            % ************************************ Stage1 分布式终止程序 ****************************************
            D=diag(diag(W));
            A=D-W;   A(A~=0) = 1;       % 邻接矩阵 (如果W不是Laplacian，则A不是邻接矩阵，那么需要把A的非零元素变成1)
            G = digraph(A');            % 注意 画图函数的矩阵元素(i,j)表示的是i到j有没有边,所以要对 A 转置一下
            dist = distances(G)';       % 返回一个矩阵dist,其元素d(i,j)表示所有节点对的最短路径距离(vj到vi的最短路径)  【后面即使W变成Wnew了，也不影响dist】
            dist(dist==inf)=0;          % 有时会出现值为inf的情况,置零即可
            [epsi_max,epsi_max_index]=max(dist,[],2);  % epsi_max是列向量,表示每个节点(dist的每行)处的偏心率;epsi_max_index表示各节点的最远距离其他节点号


            [S1_vi_ready_at_t_vec,S1_vi_stop_at_t_vec]=synchronized_termination(W,dist,epsi_max,S1_ri_max_vec);
            [S1_vi_ready_asyn_vec,S1_vi_stop_asyn_vec]=asynchronized_termination(W,dist,S1_di_vec);


        % **************************** Stage 2：行和相等的 W 矩阵转换成 列和相等的矩阵 Wnew=XWX'************************************

            X=diag(LeftVecUni_real_at_vi);     % 挑出一个节点计算出的标准左特征向量
            Wnew=X*W*inv(X);

            rng(1);
            p0_vec=ceil(rand(size(Wnew,1),1)*100-10);        % 初始向量p0_vec,随机生成一组(-10,90)之间的整数
            q0_vec(1:size(Wnew,1),1)=1;

            miu_new=sum(Wnew(:,1));          %  Wnew的列和
            p_matrix=[];   q_matrix=[];    p_matrix(:,1)=p0_vec;   q_matrix(:,1)=q0_vec;   
            t=1;      t_terminal=4*size(Wnew,1)+1;          % 在"+1"前设置迭代总次数
        %     t=1;      t_terminal=25+1;                  % 在"+1"前设置迭代总次数 (用于后续画图)
            while true
                p_matrix=[p_matrix,Wnew*p_matrix(:,t)];     % p_matrix和q_matrix用于保存p0_vec和q0_vec经过Wnew迭代的结果
                q_matrix=[q_matrix,Wnew*q_matrix(:,t)];    
                t=t+1;                
                if (t==t_terminal)       
        %             ratio_output=p_matrix./q_matrix;        % 得到p和q历史状态矩阵后,存取其ratio矩阵 (用于后续画图)
                    break;
                end
            end

        %*********************** 计算平均一致值 ******************************************
            I=eye(size(Wnew,1));
            for vi=1:size(Wnew,1)  
                eval(['eT.eT',num2str(vi),'=I(vi,:)']);
            end 

            S2_NodeFailRate=0;          %  该变量用于记录Stage 1的失败节点率
            S2_MaxTimestep=0;           %  该变量用于记录Stage 1的最大时间步
            for  vi=1:size(Wnew,1)   
                eval(['eTi=eT.eT',num2str(vi)]);
                [di,ri,average_value_at_vi,MaxTimestep] = computing_Average_consensus(Wnew,eTi,p_matrix(vi,:),q_matrix(vi,:));      % 计算平均值
                eval(['average_value_at_v',num2str(vi),'=average_value_at_vi']);     
                if vi==1
                    S2_ri_vec=ri;
                    S2_di_vec=di;
                    S2_MaxTimestep=MaxTimestep;   % 不同节点的最大时间步 取最大值
                else
                    S2_ri_vec=[S2_ri_vec;ri];       % 用向量记录Stage 2的所有的ri
                    S2_di_vec=[S2_di_vec;di];       % 用向量记录Stage 2的所有的di
                    S2_MaxTimestep=max(S2_MaxTimestep,MaxTimestep);   % 不同节点的最大时间步 取最大值
                end

                 % ****************** 测试在节点vi处,对于左特征向量x的计算是否正确 ***********
                        if miu_new==1
                            % 若W是行随机矩阵,失败包含两种情况:1.出现计算错误的节点  2. vi处的di与某个通道的rj不同(即初始条件p0属于文献中描述的零测集)
                            standard_result=roundn(sum(p0_vec)/size(Wnew,1),-2);
                            computed_result=roundn(average_value_at_vi,-2);
                            if ~isempty(find((standard_result-computed_result)~=0)) || length(find(ri~=di))>0
                                S2_NodeFailRate=S2_NodeFailRate+1;
                            end
                            if vi==size(Wnew,1)
                                S2_NodeFailRate=S2_NodeFailRate/size(Wnew,1);  % 计算节点失败率
                            end
                        else
                            % 若W不是行随机矩阵,仅统计计算出现错误的节点比例
                            standard_result=roundn(sum(p0_vec)/size(Wnew,1),-2);
                            computed_result=roundn(average_value_at_vi,-2);
                            if ~isempty(find((standard_result-computed_result)~=0))                
                                S2_NodeFailRate=S2_NodeFailRate+1;
                            end
                            if vi==size(Wnew,1)
                                S2_NodeFailRate=S2_NodeFailRate/size(Wnew,1);  % 计算节点失败率
                            end
                        end
                % *************************************************************************
 
                average_value_at_vi=[];     ri=[]; di=[];
            end

             if DatasetNumber==2
                 % 分别存储最大时间步和节点失败率.
                 % 其中,最大时间步数据用于对比: 本文的Stage 1的算法和[28]的左特征向量算法,以及本文的two stage算法和[33]的单阶段平均共识算法;
                 % 而节点失败率数据用于对比: 本文的two stage算法和[33]的单阶段平均共识算法;
                if miu_new==1
                    eval(['Output_KNNnet_Stableweight_S2_MaxTimestep.S2_MaxTimestep_of_',num2str(Num),'_nodes.MT',num2str(ni),'=S2_MaxTimestep']);
                    eval(['Output_KNNnet_Stableweight_S2_NodeFailRate.S2_NodeFailRate_of_',num2str(Num),'_nodes.NFR',num2str(ni),'=S2_NodeFailRate']);
%                     save('Output_KNNnet_Stableweight_S2_MaxTimestep','Output_KNNnet_Stableweight_S2_MaxTimestep');
                    save('Output_KNNnet_Stableweight_S2_NodeFailRate','Output_KNNnet_Stableweight_S2_NodeFailRate');  

                    S1S2_MaxTimestep=max(S1_MaxTimestep,S2_MaxTimestep);            % 两个阶段整体的最大时间步,是取Stage 1和Stage 2最大时间步的最大值
%                     
                    S1S2_NodeFailRate=max(S1_NodeFailRate,S2_NodeFailRate);            % 若阶段1的失败率是1,那么整体的失败率是1,因此不是求和 
                    eval(['Output_KNNnet_Stableweight_S1S2_MaxTimestep.S1S2_MaxTimestep_of_',num2str(Num),'_nodes.MT',num2str(ni),'=S1S2_MaxTimestep']);
                    eval(['Output_KNNnet_Stableweight_S1S2_NodeFailRate.S1S2_NodeFailRate_of_',num2str(Num),'_nodes.NFR',num2str(ni),'=S1S2_NodeFailRate']);   
%                     save('Output_KNNnet_Stableweight_S1S2_MaxTimestep','Output_KNNnet_Stableweight_S1S2_MaxTimestep');
                    save('Output_KNNnet_Stableweight_S1S2_NodeFailRate','Output_KNNnet_Stableweight_S1S2_NodeFailRate');
                else
                    eval(['Output_KNNnet_DWS1weight_S2_MaxTimestep.S2_MaxTimestep_of_',num2str(Num),'_nodes.MT',num2str(ni),'=S2_MaxTimestep']);
                    eval(['Output_KNNnet_DWS1weight_S2_NodeFailRate.S2_NodeFailRate_of_',num2str(Num),'_nodes.NFR',num2str(ni),'=S2_NodeFailRate']);
%                     save('Output_KNNnet_DWS1weight_S2_MaxTimestep','Output_KNNnet_DWS1weight_S2_MaxTimestep');
%                     save('Output_KNNnet_DWS1weight_S2_NodeFailRate','Output_KNNnet_DWS1weight_S2_NodeFailRate');  

                    S1S2_MaxTimestep=max(S1_MaxTimestep,S2_MaxTimestep);            % 两个阶段整体的最大时间步,是取Stage 1和Stage 2最大时间步的最大值
%                   
                    S1S2_NodeFailRate=max(S1_NodeFailRate,S2_NodeFailRate);            % 若阶段1的失败率是1,那么整体的失败率是1,因此不是求和 
                    eval(['Output_KNNnet_DWS1weight_S1S2_MaxTimestep.S1S2_MaxTimestep_of_',num2str(Num),'_nodes.MT',num2str(ni),'=S1S2_MaxTimestep']);
                    eval(['Output_KNNnet_DWS1weight_S1S2_NodeFailRate.S1S2_NodeFailRate_of_',num2str(Num),'_nodes.NFR',num2str(ni),'=S1S2_NodeFailRate']);   
%                     save('Output_KNNnet_DWS1weight_S1S2_MaxTimestep','Output_KNNnet_DWS1weight_S1S2_MaxTimestep');
                    save('Output_KNNnet_DWS1weight_S1S2_NodeFailRate','Output_KNNnet_DWS1weight_S1S2_NodeFailRate');           

                end
            end

            [S2_vi_ready_at_t_vec,S2_vi_stop_at_t_vec]=synchronized_termination(Wnew,dist,epsi_max,S2_ri_vec);
            [S2_vi_ready_asyn_vec,S2_vi_stop_asyn_vec]=asynchronized_termination(Wnew,dist,S2_di_vec);  

            if DatasetNumber==2
                %   存储算法的终止时间步数,该数据用于对比:本文的two stage算法(同步终止)和[33]的单阶段平均共识算法(异步终止).
                if miu_new==1
                    S2_asyn_TerminalTime=max(S2_vi_stop_asyn_vec);
                    eval(['Output_KNNnet_Stableweight_S2_asyn_TerminalTime.S2_asyn_TerminalTime_of_',num2str(Num),'_nodes.TT',num2str(ni),'=S2_asyn_TerminalTime']);
                    save('Output_KNNnet_Stableweight_S2_asyn_TerminalTime','Output_KNNnet_Stableweight_S2_asyn_TerminalTime');

                else
                    S1S2_syn_TerminalTime=max(S1_vi_stop_at_t_vec)+max(S2_vi_stop_at_t_vec);
                    eval(['Output_KNNnet_DWS1weight_S1S2_syn_TerminalTime.S1S2_syn_TerminalTime_of_',num2str(Num),'_nodes.TT',num2str(ni),'=S1S2_syn_TerminalTime']);
                    save('Output_KNNnet_DWS1weight_S1S2_syn_TerminalTime','Output_KNNnet_DWS1weight_S1S2_syn_TerminalTime');
                end
            end

            clear k ri di t t_terminal vi plot_y_vi average_value_at_vi eT eTi I D A X LeftVecUni_real_at_vi p_matrix q_matrix epsi_max_index MaxTimestep standard_result computed_result
            
        
    end
end


