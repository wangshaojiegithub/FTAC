clear;clc;
close all;
warning off;


%  **************************  (1) 首先,产生一个100节点的KNN网络,并展示效果 *********************************

%   产生规则网络
%   Num--顶点个数，Per--增加新边的概率，K--每个节点左右两边的邻居数
Num          = 100;
K            = 2;
Per          = 0;       % per=0 means that there is not random link. 
matrix_100 = func_NW_network(Num,K,Per);


% 计算相应的指标
[Cc,Cc_avg]          = func_Cluster_Coeff(matrix_100);
disp(['聚类系数为：',num2str(Cc_avg)]);
[Dds,Dds_avg,M,P_Dds]= func_Degree_Distribution(matrix_100);
disp(['平均度为：',num2str(Dds_avg)]);   
[Lens,Lens_avg]      = func_Path_Length(matrix_100);   
disp(['平均路径长度为：',num2str(Lens_avg)]); 


figure(1);  
subplot(211);   bar([1:Num],Dds);      xlabel('Node ID');  ylabel('The degree of node');   
title('Statistics on degree distribution of KNN network network with 100 nodes'); 
subplot(212);   bar([0:M],P_Dds,'r');  xlabel('The degree of node');    ylabel('Probability of node degree');
 

G = digraph(matrix_100');        %   画图,注意 画图函数的矩阵元素(i,j)表示的是i到j有没有边,所以要对 A 转置一下
figure(2);  
p_graph=plot(G,'Layout','circle','EdgeColor',[0 0.4 0.8],'LineWidth',0.5,'NodeColor','k','NodeFontSize',8);   % 画出权值图 

[average_loopslength,MaximumLoopLength] = func_averageLML(G,matrix_100,p_graph);             % 计算平均LML，并在计算过程中用红色边标出最大的loop
disp(['100节点的规则网络的LML长度为：',num2str(MaximumLoopLength)]); 
disp(['100节点的规则网络的LMLs平均长度为：',num2str(average_loopslength)]); 



% % **********************  (2)  网络规模设定：50,100,150,...,950,1000.每个网络规模下生成20个网络以观察平均表现 *********************
for n=50:50:1000
    Num          = n;
    for ni=1:20
        Adja_matrix  = func_NW_network(Num,K,Per);
%         while true
            Adja_matrix_temp=Adja_matrix;
            for i=1:Num 
                for j=(i+1):Num
                    if Adja_matrix_temp(i,j)~=0 && rand()<0.1
                       Adja_matrix_temp(i,j)=0;               % 以概率0.1删除原图邻接矩阵中上三角的一些边
                    end
                end
            end
%             % 确保随机删除连边后的图仍然是强连通的
%                        Gtemp = digraph(Adja_matrix_temp');        
%                        dist = distances(Gtemp)';            
%             L=diag(sum(Adja_matrix_temp,2))-Adja_matrix_temp; eigenvalues=eig(L);  
%             eigenvalues_003=roundn(eigenvalues,-3);   
%             if length(find(eigenvalues_003==0))==1 && isempty(dist(dist==inf))
%                 dist=[]; Gtemp=[];
%                 Adja_matrix=Adja_matrix_temp;
%                 break;
%             end
%         end
        eval(['KNN_matrices.KNN_matrices_of_',num2str(Num),'_nodes.A',num2str(ni),'=Adja_matrix']);       
    end    
end
save ('KNN_matrices','KNN_matrices');     % 所有数据保存到一个结构体


% ************************ (3) 导入保存的邻接矩阵数据 BA_matrices.mat, 完成不同规模网络的相关统计 ****************************
load('KNN_matrices.mat');
% load('KNNnet_AveLengthofLMLs.mat');

for n=50:50:1000
    Num          = n;
    for ni=1:20
        eval(['matrix_temp=KNN_matrices.KNN_matrices_of_',num2str(Num),'_nodes.A',num2str(ni)]); % 相当于matrix_temp=BA_matrices.BA_matrices_of_200_nodes.A2;
        
        G_temp = digraph(matrix_temp');    figure(3);  
        p_graph_temp=plot(G_temp,'Layout','force','EdgeColor',[0 0.4 0.8],'LineWidth',0.5,'NodeColor','k','NodeFontSize',8);   
        [AverageLength_of_LMLs,MaxLoopLength] = func_averageLML(G_temp,matrix_temp,p_graph_temp);      % 计算平均LML,并在计算过程中用红色边标出最大的loop
        
        eval(['KNNnet_AveLengthofLMLs.KNNnet_AveLengthofLMLs_',num2str(Num),'_nodes.A',num2str(ni),'=AverageLength_of_LMLs']);
%         eval(['BAnet_MaxLoopLength.BAnet_MaxLoopLength_',num2str(Num),'_nodes.A',num2str(ni),'=MaxLoopLength']);
        
        save ('KNNnet_AveLengthofLMLs','KNNnet_AveLengthofLMLs');
%         save ('BAnet_MaxLoopLength','BAnet_MaxLoopLength');
    end
end



