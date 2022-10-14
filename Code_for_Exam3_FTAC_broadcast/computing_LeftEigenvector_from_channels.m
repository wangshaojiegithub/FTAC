function [di_real,ri_max,ri_max_source_vec,left_vi_real_Uni,left_vi_Uni]=computing_LeftEigenvector_from_channels(vi,ei,miu,W,p1_vec,p2_vec,p3_vec,p4_vec,q0_vec)

t=1;          t_terminal=30+1;                   % 在"+1"前设置迭代总次数
p_matrix=[];  q_matrix=[]; O_alpha_matrix=[]; O_delta_matrix=[];

%********** 开始计算：在vi处计算标准的di, delta_i向量和left vec(表示成left_vi_real_Uni)**********
while_flag=1;    t=1;
O_alpha_matrix(1,:)=ei;
O_delta_matrix(1,:)=ei*(W-miu*eye(size(W,1)));
while(while_flag==1)
    O_alpha_matrix=[O_alpha_matrix;O_alpha_matrix(t,:)*W];
    O_delta_matrix=[O_delta_matrix;O_delta_matrix(t,:)*W];
    t=t+1;             % 表示已经具有t个连续数据(向量),其实迭代次数(时间)是t-1    
    if (t==t_terminal)
        while_flag=0; 
    end
end
% 标准delta_i向量
while_stop_flag=0;
AAA_temp=[];
i=2;
while (while_stop_flag==0 & i<=size(O_delta_matrix,1))
    AAA_temp=O_delta_matrix(1:i,:);
    if (rank(AAA_temp)<size(AAA_temp,1))
        delta_i_real=null(AAA_temp','r')';    % 可用于计算标准的delta_i 
        di_real=rank(AAA_temp');              % 标准的di,直接作为输出数据
        while_stop_flag=1;
    else
        AAA_temp=[];
        i=i+1;
    end        
end
% 在vi处计算出的标准化的左特征向量(与初值无关)
left_vi_real=delta_i_real*O_alpha_matrix(1:size(delta_i_real,2),:); 
left_vi_real_Uni=left_vi_real/sum(left_vi_real);    % 直接作为输出数据
%     eval(['left_',num2str(vi),'_real_Uni=left_vi_real_Uni']);
%************************************************************************

% 在vi处,接下来的每个ci表示第i个channel(共n个)
for  cj=1:size(W,1)    
    eval(['pi_vec=p',num2str(cj),'_vec']);
    t=1;   t_terminal=30+1;
    p_matrix(:,1)=pi_vec;
    q_matrix(:,1)=q0_vec;
    while_flag=1;
    while(while_flag==1)
        p_matrix=[p_matrix,W*p_matrix(:,t)];
        q_matrix=[q_matrix,W*q_matrix(:,t)];
        t=t+1;             % 表示已经具有t个连续数据(向量),其实迭代次数(时间)是t-1 
        if (t==t_terminal)
            while_flag=0; 
        end
    end
    pi_series=ei*p_matrix;
    qi_series=ei*q_matrix;

    pi_series_dif=[];
    for i=1:size(pi_series,2)-1
        pi_series_dif(1,i)=pi_series(1,i+1)-miu*pi_series(1,i);%由p矩阵第i行构建一个差分向量
    end    
    pi_delta_FDH=[];  % 差分Hankel矩阵
    while_stop_flag=0;
    ri_j=1;         % 用来记录vi处计算第j通道值的 差分FDH矩阵的秩

    k=8;            %   差分向量构建的Hankel矩阵要从k=2维(或di或n维)开始考察
    while (while_stop_flag==0 & 2*k-1<=size(pi_series_dif,2))
        for i=1:k
            pi_delta_FDH(i,:)=pi_series_dif(1,i:k+i-1);
        end 
        if (rank(pi_delta_FDH)<k)
            while_stop_flag=1;  % 准备跳出while循环
        else
            pi_delta_FDH=[];    % 注意如果该Hankel不是缺陷的，要在下一次循环清空，以免矩阵维数被固定
            k=k+1;
        end   
    end    
    
    ri_j=rank(pi_delta_FDH);   % 其实现在的Hankel矩阵的秩是r=k-1。
    delta_i_j=null(pi_delta_FDH(1:ri_j+1,1:ri_j+1),'r')';
    left_i_j=(delta_i_j*pi_series(1:size(delta_i_j,2))')/(delta_i_j*qi_series(1:size(delta_i_j,2))'); % 标准的(归一化)左特征向量第i个分量
    
    %   记录要输出的数据
    if cj==1
        ri_max=ri_j;
        ri_max_source_vec=ri_j;
        
        left_vi_Uni=left_i_j;
    else
        ri_max=max([ri_max,ri_j]);
        ri_max_source_vec=[ri_max_source_vec,ri_j];
        
        left_vi_Uni=[left_vi_Uni,left_i_j];
    end
      
    p_matrix=[];  q_matrix=[]; O_alpha_matrix=[]; O_delta_matrix=[];     delta_i_j=[];
end 
end

