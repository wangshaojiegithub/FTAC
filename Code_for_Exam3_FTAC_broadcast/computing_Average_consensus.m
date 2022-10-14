function [ri,average_value_at_vi,delta_i_vec,Qi]=computing_Average_consensus(pi_series,qi_series,miu_new)
    pi_series_dif=[]; qi_series_dif=[];
    pi_delta_FDH=[];   qi_delta_FDH=[]; 
    for i=1:size(pi_series,2)-1
        pi_series_dif(1,i)=pi_series(1,i+1)-miu_new*pi_series(1,i);     % ��p�����i�й���һ���������
    end   
    for i=1:size(qi_series,2)-1
        qi_series_dif(1,i)=qi_series(1,i+1)-miu_new*qi_series(1,i);%��p�����i�й���һ���������
    end  
   
    % ���Hankel����   
    while_stop_flag=0;
    ri=1;               % ������¼First Hankel�������    

    k=4;                % �������������Hankel����Ҫ��k=2ά(��di��nά)��ʼ���� ��ע�⣺k=4(������8Ŷ)��
    while (while_stop_flag==0 & 2*k-1<=size(pi_series_dif,2))
        for i=1:k
            pi_delta_FDH(i,:)=pi_series_dif(1,i:k+i-1);
            qi_delta_FDH(i,:)=qi_series_dif(1,i:k+i-1);
        end 
        
    %********** ��1����򣬻�����׵����޷����ƽ��ֵ����ζ��p0����pi_series���� ������⼯***************
%         if (rank(pi_delta_FDH)<k)
%             while_stop_flag=1;  % ׼������whileѭ��
%             ri=rank(pi_delta_FDH);  
%             
% %             delta_i_vec=null(pi_delta_FDH(1:ri+1,1:ri+1),'r')';
%             
%             delta_i_vec_whole=null(pi_delta_FDH,'r')';
%             delta_i_vec=delta_i_vec_whole(1,1:ri+1);
%             
%             average_value_at_vi=(delta_i_vec*pi_series(1:size(delta_i_vec,2))')/(delta_i_vec*qi_series(1:size(delta_i_vec,2))'); 
%         else
%             pi_delta_FDH=[];    % ע�������Hankel����ȱ�ݵģ�Ҫ����һ��ѭ����գ��������ά�����̶�
%             k=k+1;
%         end  
     %****************************************************************************************************
     %********** ��2�����ͬʱ�ж�p��q������ֵ���Hankel���ȣ�ȡ��������ڼ��㣬�ɼ����޷����ƽ��ֵ��״�� *** 
        if (rank(pi_delta_FDH)<k & rank(qi_delta_FDH)<k)                       
            if rank(pi_delta_FDH)>=rank(qi_delta_FDH)
                ri=rank(pi_delta_FDH);             
                delta_i_vec=null(pi_delta_FDH(1:ri+1,1:ri+1),'r')';       % !!!ע��:д��pi_delta_FDH(1:ri+1,1:ri+1)������pi_delta_FDH
            else
                ri=rank(qi_delta_FDH);             
                delta_i_vec=null(qi_delta_FDH(1:ri+1,1:ri+1),'r')';
            end          
            % ����ж�:q_matrix_vi_row*delta_vec_temp �Ƿ�Ϊ0������ǣ���������Hankelά�����Ա���ratio��ĸΪ��.(��ʵҲֻ������ǰ���ȣ�������ȫ�������)
%             if qi_series(1:size(delta_i_vec,1))*delta_i_vec(:,1)==0
%                 pi_delta_FDH=[]; qi_delta_FDH=[];   delta_i_vec=[];         
%                 k=k+1; 
%             else
                while_stop_flag=1;         % ׼������whileѭ��
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