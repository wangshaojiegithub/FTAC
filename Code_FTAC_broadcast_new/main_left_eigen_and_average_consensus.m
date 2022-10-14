clear;clc;
% **************Stage 1：计算行和相等的 W=L 的左特征向量************
W=[  2 0  0  0  0  0  0 -1;
    -1 1 -1  0  0  0  0  0;
    0  0  3  0  0  0  0 -1;
    0  0 -1  1 -1 0  0  0;
    0  0  0  0  3  0  0 -1;
    0  0  0  0 -1  1 -1  0;
    0  0  0  0  0  0  1 -1;
   -1 -1 -1 -1 -1 -1 0  4];      % 列和为0  （例1中Algorithm 2 用的矩阵）

% W=[ 1/2  0     0   0    0    0    0   1/8;
%     1/4  1/2  1/6  0    0    0    0    0;
%     0    0    1/2  0    0    0    0   1/8;
%     0    0    1/6  1/2  1/6  0    0    0;
%     0    0    0    0    1/2  0    0   1/8;
%     0    0    0    0    1/6  1/2  1/2  0;
%     0    0    0    0    0    0    1/2 1/8;
%     1/4 1/2   1/6  1/2  1/6  1/2  0   1/2];      % 列和为1 （修改稿中例1复现[28]用的矩阵，其实结果和之前一样，因为两个矩阵其实都是L型的）

    W=W';                           % 转成行和相等的矩阵
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
                    standard_result=roundn(LeftVecUni_real_at_vi,-4);
                    computed_result=roundn(LeftVecUni_at_vi,-4);
                    if ~isempty(find((standard_result-computed_result)~=0)) || length(find(ri_max_source_vec~=di_real))>0
                        S1_NodeFailRate=S1_NodeFailRate+1;
                    end
                    if vi==size(W,1)
                        S1_NodeFailRate=S1_NodeFailRate/size(W,1); % 计算节点失败率
                    end
                else
                    % 若W不是行随机矩阵,仅统计计算出现错误的节点比例
                    standard_result=roundn(LeftVecUni_real_at_vi,-4);
                    computed_result=roundn(LeftVecUni_at_vi,-4);
                    if ~isempty(find((standard_result-computed_result)~=0))                
                        S1_NodeFailRate=S1_NodeFailRate+1;
                    end
                    if vi==size(W,1)
                        S1_NodeFailRate=S1_NodeFailRate/size(W,1); % 计算节点失败率
                    end
                end
            % *************************************************************************
        eval(['S1_ri_max_vec_source.r',num2str(vi),'max_source=ri_max_source_vec']); % 用于检查ri_max的每个元素的来源,保存为一个结构体

        eval(['S1_LeftVecUni_real_at_v',num2str(vi),'=LeftVecUni_real_at_vi']);
        eval(['S1_LeftVecUni_at_v',num2str(vi),'=LeftVecUni_at_vi']);      
    end
    save('S1_MaxTimestep','S1_MaxTimestep');
    save('S1_NodeFailRate','S1_NodeFailRate');
    clearvars -except W LeftVecUni_real_at_vi S1_MaxTimestep S1_FailRate_of_nodes S1_di_vec S1_ri_max_vec S1_ri_max_vec_source

  % ************************************ Stage1 分布式终止程序 *********************************************************
    D=diag(diag(W));
    A=D-W;   A(A~=0) = 1;       % 邻接矩阵 (如果W不是Laplacian，则A不是邻接矩阵，那么需要把A的非零元素变成1)
    G = digraph(A');            % 注意 画图函数的矩阵元素(i,j)表示的是i到j有没有边,所以要对 A 转置一下
    dist = distances(G)';       % 返回一个矩阵dist,其元素d(i,j)表示所有节点对的最短路径距离(vj到vi的最短路径)  【后面即使W变成Wnew了，也不影响dist】
    [epsi_max,epsi_max_index]=max(dist,[],2);  % epsi_max是列向量,表示每个节点(dist的每行)处的偏心率;epsi_max_index表示各节点的最远距离其他节点号

    figure(1);  p_graph=plot(G,'EdgeColor',[0 0.0 0.0],'LineWidth',1,'NodeColor','k','NodeFontSize',12);        % 画出权值图
    set(gcf,'color','none');

    [S1_vi_ready_at_t_vec,S1_vi_stop_at_t_vec]=synchronized_termination(W,dist,epsi_max,S1_ri_max_vec);
    [S1_vi_ready_asyn_vec,S1_vi_stop_asyn_vec]=asynchronized_termination(W,dist,S1_di_vec);

    figure(2);                  % ************* 画出终止Stage 1迭代的终止时间图
    plot_y_vi=1:1:size(W,1);  
    syn_stop_S1_plot =  plot(S1_vi_ready_at_t_vec,plot_y_vi,'.','color','m','MarkerSize',10); hold on;
    syn_stop_S1_plot =  plot(S1_vi_stop_at_t_vec,plot_y_vi,'*','color','m','MarkerSize',5); hold on;
    asyn_stop_S1_plot = plot(S1_vi_ready_asyn_vec,plot_y_vi,'+',S1_vi_stop_asyn_vec,plot_y_vi,'*','color','b','MarkerSize',4);hold off;

    xlabel('Time step k (iteration)');ylabel('Agent ID (i)'); title('Time step for terminating iterations in Stage 1');
    % axis([0 25 0 9]); set(gca,'XTick',[0:2:25]); set(gca,'YTick',[0:1:9]);
    axis([0 2*size(W,1) 0 size(W,1)+1]); set(gca,'XTick',[0:2:2*size(W,1)]); set(gca,'YTick',[0:1:size(W,1)+1]);
    legend('Ready (Algorithm 2)','Stop    (Algorithm 2)','Ready (algorithm in [33])','Stop    (algorithm in [33])');
  % **************************************************************************************************************

% **************************** Stage 2：行和相等的 W 矩阵转换成 列和相等的矩阵 Wnew=XWX'************************************

    X=diag(LeftVecUni_real_at_vi);     % 挑出一个节点计算出的标准左特征向量
    Wnew=X*W*inv(X);

% ***********   这段代码用于测试:例1中Stage 2的对比算法[27],取消注释即可运行     **********
    %               虽然本程序能计算最终结果,但一些节点的di≠ri,这在[27]的算法下属于计算错误.
    
%     clearvars -except dist epsi_max S2_ri_vec S1_di_vec
%     Wnew=[ 1/3  0     0   0    0    0    0   1/5;
%         1/3  1/2  1/4  0    0    0    0    0;
%         0    0    1/4  0    0    0    0   1/5;
%         0    0    1/4  1/2  1/4  0    0    0;
%         0    0    0    0    1/4  0    0   1/5;
%         0    0    0    0    1/4  1/2  1/2  0;
%         0    0    0    0    0    0    1/2 1/5;
%         1/3 1/2   1/4  1/2  1/4  1/2  0   1/5];
    
% *************************************************************************************

    % p0_vec=[1;0;4;0;-9;0;0;3];            % 属于零测集,可以通过"同时loss rank"条件来减弱这种状况(初稿的例2)
    p0_vec=[1;4;4;3;-9;2;5;3];              % 不属于零测集 (初稿的例1)
    q0_vec([1:size(Wnew,1)],1)=1;
    
    miu_new=sum(Wnew(:,1));          %  Wnew的列和
    p_matrix=[];   q_matrix=[];    p_matrix(:,1)=p0_vec;   q_matrix(:,1)=q0_vec;   
    % t=1;      t_terminal=3*size(Wnew,1)+1;          % 在"+1"前设置迭代总次数
    t=1;      t_terminal=25+1;          % 在"+1"前设置迭代总次数 (用于后续画图)
    while true
        p_matrix=[p_matrix,Wnew*p_matrix(:,t)];     % p_matrix和q_matrix用于保存p0_vec和q0_vec经过Wnew迭代的结果
        q_matrix=[q_matrix,Wnew*q_matrix(:,t)];    
        t=t+1;                
        if (t==t_terminal)       
            ratio_output=p_matrix./q_matrix;        % 得到p和q历史状态矩阵后,存取其ratio矩阵 (用于后续画图)
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
                    standard_result=roundn(sum(p0_vec)/size(Wnew,1),-4);
                    computed_result=roundn(average_value_at_vi,-4);
                    if ~isempty(find((standard_result-computed_result)~=0)) || length(find(ri~=di))>0
                        S2_NodeFailRate=S2_NodeFailRate+1;
                    end
                    if vi==size(Wnew,1)
                        S2_NodeFailRate=S2_NodeFailRate/size(Wnew,1);  % 计算节点失败率
                    end
                else
                    % 若W不是行随机矩阵,仅统计计算出现错误的节点比例
                    standard_result=roundn(sum(p0_vec)/size(Wnew,1),-4);
                    computed_result=roundn(average_value_at_vi,-4);
                    if ~isempty(find((standard_result-computed_result)~=0))                
                        S2_NodeFailRate=S2_NodeFailRate+1;
                    end
                    if vi==size(Wnew,1)
                        S2_NodeFailRate=S2_NodeFailRate/size(Wnew,1);  % 计算节点失败率
                    end
                end
        % *************************************************************************
                
        for k=1:2*(ri+1)-1
            ratio_output(vi,k)=nan;                     % 先将ratio_output的前2*(ri+1)个值设置为NAN,用于后续画图
        end
        for k=2*(ri+1):size(ratio_output(vi,:),2)
            ratio_output(vi,k)=average_value_at_vi;     % 再将2*(ri+1)及其之后的ratio_output值设置为average value
        end     
        average_value_at_vi=[];     ri=[]; di=[];
    end
    save('S2_MaxTimestep','S2_MaxTimestep');
    save('S2_NodeFailRate','S2_NodeFailRate');

    [S2_vi_ready_at_t_vec,S2_vi_stop_at_t_vec]=synchronized_termination(Wnew,dist,epsi_max,S2_ri_vec);
    [S2_vi_ready_asyn_vec,S2_vi_stop_asyn_vec]=asynchronized_termination(Wnew,dist,S2_di_vec);  

    for  vi=1:size(Wnew,1)   
        for k=S2_vi_stop_at_t_vec+2:size(ratio_output(vi,:),2)
            % 要展示同步终止,将 S2_vi_stop_at_t_vec+2 及其之后的p,q和ratio_output值设置为NAN
            p_matrix(vi,k)=nan;
            q_matrix(vi,k)=nan;
            ratio_output(vi,k)=nan;  
        end
    end

%*************** 画出每个节点的两个迭代(p和q)的历史状态图以及ratio图 ************
    color = colormap(jet(size(Wnew,1))); 
    plot_t=0:1:t_terminal-1;                % 绘制历史曲线（共t_terminal个时刻,包含0时刻）
    for vi=1:size(Wnew,1)
        figure(3);
        subplot(3,1,1);
        p_plot=plot(plot_t,p_matrix(vi,1:t_terminal),'-','color',color(vi,:)); % 绘制第一个迭代的历史状态.例如plot_t=100，共100个历史时刻（99次迭代）
        ylabel('p-value');
        xlim([0 t_terminal-1]);
        title('The "ratio iterations" in Stage 2');
        hold on;grid on;

        subplot(3,1,2)
        q_plot=plot(plot_t,q_matrix(vi,1:t_terminal),'-','color',color(vi,:)); % 绘制第二个迭代的历史状态
        ylabel('q-value');
        xlim([0 t_terminal-1]);
        hold on;grid on;

        subplot(3,1,3)
        r_plot=plot(plot_t,ratio_output(vi,1:t_terminal),'-','color',color(vi,:));
        xlabel('Time step k (iteration)');ylabel('Ratio value');
    %     axis([0 t_terminal-1 -0.13 -0.12]);        % 限定一下横坐标和纵坐标范围（用于原稿的例2）
            axis([0 t_terminal-1 1.60 1.65]);        % 限定一下横坐标和纵坐标范围（用于原稿的例1）
        hold on;grid on;
    end
    legend('{\it{v}}1','{\it{v}}2','{\it{v}}3','{\it{v}}4','{\it{v}}5','{\it{v}}6','{\it{v}}7','{\it{v}}8');

clear k ri di t t_terminal vi while_flag plot_y_vi average_value_at_vi eT eTi I D A p_matrix q_matrix epsi_max_index MaxTimestep standard_result computed_result
