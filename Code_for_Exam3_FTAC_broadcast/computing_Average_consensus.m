function [ri,average_value_at_vi,delta_i_vec,Qi]=computing_Average_consensus(pi_series,qi_series,miu_new)
    pi_series_dif=[]; qi_series_dif=[];
    pi_delta_FDH=[];   qi_delta_FDH=[]; 
    for i=1:size(pi_series,2)-1
        pi_series_dif(1,i)=pi_series(1,i+1)-miu_new*pi_series(1,i);     % 由p矩阵第i行构建一个差分向量
    end   
    for i=1:size(qi_series,2)-1
        qi_series_dif(1,i)=qi_series(1,i+1)-miu_new*qi_series(1,i);%由p矩阵第i行构建一个差分向量
    end  
   
    % 差分Hankel矩阵   
    while_stop_flag=0;
    ri=1;               % 用来记录First Hankel矩阵的秩    

    k=4;                % 差分向量构建的Hankel矩阵要从k=2维(或di或n维)开始考察 【注意：k=4(而不是8哦)】
    while (while_stop_flag==0 & 2*k-1<=size(pi_series_dif,2))
        for i=1:k
            pi_delta_FDH(i,:)=pi_series_dif(1,i:k+i-1);
            qi_delta_FDH(i,:)=qi_series_dif(1,i:k+i-1);
        end 
        
    %********** 第1版程序，会更容易导致无法算出平均值，意味着p0对于pi_series而言 属于零测集***************
%         if (rank(pi_delta_FDH)<k)
%             while_stop_flag=1;  % 准备跳出while循环
%             ri=rank(pi_delta_FDH);  
%             
% %             delta_i_vec=null(pi_delta_FDH(1:ri+1,1:ri+1),'r')';
%             
%             delta_i_vec_whole=null(pi_delta_FDH,'r')';
%             delta_i_vec=delta_i_vec_whole(1,1:ri+1);
%             
%             average_value_at_vi=(delta_i_vec*pi_series(1:size(delta_i_vec,2))')/(delta_i_vec*qi_series(1:size(delta_i_vec,2))'); 
%         else
%             pi_delta_FDH=[];    % 注意如果该Hankel不是缺陷的，要在下一次循环清空，以免矩阵维数被固定
%             k=k+1;
%         end  
     %****************************************************************************************************
     %********** 第2版程序，同时判断p和q两个差分迭代Hankel的秩，取最大者用于计算，可减轻无法算出平均值的状况 *** 
        if (rank(pi_delta_FDH)<k & rank(qi_delta_FDH)<k)                       
            if rank(pi_delta_FDH)>=rank(qi_delta_FDH)
                ri=rank(pi_delta_FDH);             
                delta_i_vec=null(pi_delta_FDH(1:ri+1,1:ri+1),'r')';       % !!!注意:写成pi_delta_FDH(1:ri+1,1:ri+1)而不是pi_delta_FDH
            else
                ri=rank(qi_delta_FDH);             
                delta_i_vec=null(qi_delta_FDH(1:ri+1,1:ri+1),'r')';
            end          
            % 最后判断:q_matrix_vi_row*delta_vec_temp 是否为0。如果是，则再增加Hankel维数，以避免ratio分母为零.(其实也只缓解提前降秩，不能完全解决问题)
%             if qi_series(1:size(delta_i_vec,1))*delta_i_vec(:,1)==0
%                 pi_delta_FDH=[]; qi_delta_FDH=[];   delta_i_vec=[];         
%                 k=k+1; 
%             else
                while_stop_flag=1;         % 准备跳出while循环
%             end
        else
            pi_delta_FDH=[];      qi_delta_FDH=[];            
            k=k+1; 
        end  
    %***************************************************************************************************   
    end    
    average_value_at_vi=(delta_i_vec*pi_series(1:size(delta_i_vec,2))')/(delta_i_vec*qi_series(1:size(delta_i_vec,2))');   
    
    Qi=delta_i_vec*qi_series(1:size(delta_i_vec,2))';
end