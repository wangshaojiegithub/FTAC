clear;clc;
% **************Stage 1：计算行和相等的 W=L 的左特征向量************
W2=[1/2 0 1/2  0;
    0   1  0   0;
    0   0 1/2 1/2;
    0   0  0   1];
W1=[1/2 0  0  1/2;
    0   1  0   0;
    0   0  1   0;
    0  1/2 0  1/2];
W0=[1   0   0   0;
    1/3 1/3 1/3 0;
    0   0   1   0;
    0   0   0   1];
W=W2*W1*W0;         % 积图，列和为1

p1_vec=[1;0;0;0];
p2_vec=[0;1;0;0];
p3_vec=[0;0;1;0];
p4_vec=[0;0;0;1];

q0_vec=[1;1;1;1];
e1=[1;0;0;0]';
e2=[0;1;0;0]';
e3=[0;0;1;0]';
e4=[0;0;0;1]';

miu=sum(W(1,:));                    %求行和
% miu=sum(W(:,1));                    %求列和
for  vi=1:size(W,1)    
    eval(['ei=e',num2str(vi)]);
    [di_real,ri_max,ri_max_source_vec,left_vi_real_Uni,left_vi_Uni]=computing_LeftEigenvector_from_channels(vi,ei,miu,W,p1_vec,p2_vec,p3_vec,p4_vec,q0_vec);
%     eval(['S1_d',num2str(vi),'=di_real']);
%     eval(['S1_r',num2str(vi),'_max=ri_max']);
    % 以列向量形式表达:S1_di_vec和S1_ri_max_vec
    if vi==1
        S1_di_vec=di_real;
        S1_ri_max_vec=ri_max;
    else
        S1_di_vec=[S1_di_vec;di_real];
        S1_ri_max_vec=[S1_ri_max_vec;ri_max];
    end    
    eval(['S1_r',num2str(vi),'max_source_vec=ri_max_source_vec']); % 用于检查ri_max的每个元素的来源
    
    eval(['S1_left_v',num2str(vi),'_real_Uni=left_vi_real_Uni']);
    eval(['S1_left_v',num2str(vi),'_Uni=left_vi_Uni']);      
end
clear ei di_real ri_max vi left_vi_real_Uni left_vi_Uni
clear p1_vec p2_vec p3_vec p4_vec 
clear e1 e2 e3 e4 
clear S1_left_v2_real_Uni S1_left_v2_Uni S1_left_v3_real_Uni S1_left_v3_Uni S1_left_v4_real_Uni S1_left_v4_Uni 

% ****************** 分布式终止程序 ********************************
D=diag(diag(W));
A=D-W;   A(A~=0) = 1;       % 邻接矩阵 (如果W不是Laplacian，则A不是邻接矩阵，那么需要把A的非零元素变成1)
G = digraph(A');            % 画图,注意 画图函数的矩阵元素(i,j)表示的是i到j有没有边,所以要对 A 转置一下
dist = distances(G)';       % 返回一个矩阵dist,其元素d(i,j)表示所有节点对的最短路径距离(vj到vi的最短路径)  【后面即使W变成Wnew了，也不影响dist】
[epsi_max,epsi_max_index]=max(dist,[],2);  % epsi_max是一个列向量,表示每个节点(dist的每行)处的偏心率;epsi_max_index表示各节点的最远距离其他节点号

figure(1);
% p_graph=plot(G,'EdgeLabel',G.Edges.Weight);                   % 画出权值图
p_graph=plot(G,'EdgeColor',[0 0.0 0.0],'LineWidth',1,'NodeColor','k','NodeFontSize',12);                   % 画出权值图
set(gcf,'color','none');
% set(gca,'color','none');
% [path_max_3,distance_of_path_max_3] = shortestpath(G,epsi_max_index(3),3);  % 画到第3个节点的最长路径
% highlight(p_graph,path_max_3,'EdgeColor','r');                % 高亮第3个节点的最长路径

[S1_vi_ready_at_t_vec,S1_vi_stop_at_t_vec]=synchronized_termination(W,dist,epsi_max,S1_ri_max_vec);
[S1_vi_ready_asyn_vec,S1_vi_stop_asyn_vec]=asynchronized_termination(W,dist,S1_di_vec);

figure(2);                  % ************* 画出终止Stage 1迭代的终止时间图
plot_y_vi=1:1:size(W,1);  
syn_stop_S1_plot =  plot(S1_vi_ready_at_t_vec,plot_y_vi,'.','color','m','MarkerSize',10); hold on;
syn_stop_S1_plot =  plot(S1_vi_stop_at_t_vec,plot_y_vi,'*','color','m','MarkerSize',5); hold on;
asyn_stop_S1_plot = plot(S1_vi_ready_asyn_vec,plot_y_vi,'+',S1_vi_stop_asyn_vec,plot_y_vi,'*','color','b','MarkerSize',4);hold off;

xlabel('Time step k (iteration)');ylabel('Agent ID (i)'); title('Time step for terminating iterations in Stage 1');
axis([0 25 0 5]); set(gca,'XTick',[0:2:25]); set(gca,'YTick',[0:1:5]);
legend('Ready (Algorithm 2)','Stop    (Algorithm 2)','Ready (algorithm in [33])','Stop    (algorithm in [33])');
% **********************************************************************

% **************************** Stage 2：行和相等的 W 矩阵转换成 列和相等的矩阵 Wnew=XWX'************************************
X=diag(S1_left_v1_real_Uni); % 挑出一个节点计算出的标准左特征向量
Wnew=X*W*inv(X);

miu_new=sum(Wnew(:,1));                % 求Wnew的列和

p0_vec=[5;-7;1;3];

p_matrix=[];   q_matrix=[];
p_matrix(:,1)=p0_vec;   
q_matrix(:,1)=q0_vec;   

t=1;      t_terminal=25+1;          % 【在"+1"前设置迭代总次数】
while_flag=1;       % 矩阵p_matrix和q_matrix存储p0_vec和q0_vec经过Wnew的迭代结果
while (while_flag==1)
    p_matrix=[p_matrix,Wnew*p_matrix(:,t)];
    q_matrix=[q_matrix,Wnew*q_matrix(:,t)];    
    t=t+1;                
    if (t==t_terminal)
        while_flag=0; 
        ratio_output=p_matrix./q_matrix; % 得到p和q历史状态矩阵后,存取其ratio矩阵
    end
end

%*********************** 计算平均一致值 ******************************************
for  vi=1:size(Wnew,1)   
    [ri,average_value_at_vi,delta_i_vec,Qi] = computing_Average_consensus(p_matrix(vi,:),q_matrix(vi,:),miu_new);      
    eval(['average_value_at_v',num2str(vi),'=average_value_at_vi']);    
    
    eval(['delta_v',num2str(vi),'_vec=delta_i_vec']); % 记录delta_i_vec用于噪声网络先验知识
    eval(['Q',num2str(vi),'=Qi']);                   % 记录Qi用于噪声网络先验知识
    
    if vi==1
        S2_ri_vec=ri;
    else
        S2_ri_vec=[S2_ri_vec;ri];       % 用向量记录Stage 2的所有的ri
    end
    
    for k=1:2*(ri+1)-1
        ratio_output(vi,k)=nan;         % 先将ratio_output的前2*(ri+1)个值设置为NAN
    end
    for k=2*(ri+1):size(ratio_output(vi,:),2)
        ratio_output(vi,k)=average_value_at_vi;  % 再将2*(ri+1)及其之后的ratio_output值设置为average value
    end 
    
    average_value_at_vi=[];     ri=[];
end

[S2_vi_ready_at_t_vec,S2_vi_stop_at_t_vec]=synchronized_termination(Wnew,dist,epsi_max,S2_ri_vec);
[S2_vi_ready_asyn_vec,S2_vi_stop_asyn_vec]=asynchronized_termination(Wnew,dist,S1_di_vec);  % 注意:S1_di_vec等于"S2_di_vec"
for  vi=1:size(Wnew,1)   
    for k=S2_vi_stop_at_t_vec+2:size(ratio_output(vi,:),2)
        % 因为要展示同步终止，所以将 S2_vi_stop_at_t_vec+2 及其之后的p, q和ratio_output值设置为NAN
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

        axis([0 t_terminal-1 0.3 0.7]);        % 限定一下横坐标和纵坐标范围
    hold on;grid on;
    end
legend('{\it{v}}1','{\it{v}}2','{\it{v}}3','{\it{v}}4');

clear k ri t t_terminal vi while_flag plot_y_vi average_value_at_vi
