function [di_real,ri,average_value_at_vi,maxtimestep]=computing_Average_consensus(Wnew,eTi,pi_series,qi_series)

%********************** 在vi处计算标准的di,delta_i向量 **********************
    O_delta_matrix=[];    maxtimestep_vec=[];
    miu_new=sum(Wnew(:,1));                                % 求列和
    t=1;    t_terminal=2*size(Wnew,1)+1;                   % 在"+1"前设置迭代总次数
    O_delta_matrix(1,:)=eTi*(Wnew-miu_new*eye(size(Wnew,1)));
    while true   
        O_delta_matrix=[O_delta_matrix;O_delta_matrix(t,:)*Wnew];
        t=t+1;             % 表示已经具有t个连续数据(向量),其实迭代次数(时间)是t-1    
        if (t==t_terminal)
            break; 
        end
    end
    AAA_temp=[];
    i=2;
    while i<=size(O_delta_matrix,1)
        AAA_temp=O_delta_matrix(1:i,:);
        if (rank(AAA_temp)<size(AAA_temp,1))
    %         delta_i_real=null(AAA_temp','r')';    % 可用于计算标准的delta_i 
            di_real=rank(AAA_temp');              % 标准的di,直接作为输出数据
            break;
        else
            AAA_temp=[];
            i=i+1;
        end        
    end 
% **************************************************************************

    pi_series_dif=[]; qi_series_dif=[];
    pi_delta_FDH=[];   qi_delta_FDH=[]; 
    for i=1:size(pi_series,2)-1
        pi_series_dif(1,i)=pi_series(1,i+1)-miu_new*pi_series(1,i);     % 由p矩阵第i行构建一个差分向量
    end   
    for i=1:size(qi_series,2)-1
        qi_series_dif(1,i)=qi_series(1,i+1)-miu_new*qi_series(1,i);     % 由p矩阵第i行构建一个差分向量
    end  
   
    ri=1;               % 用来记录First Hankel矩阵的秩    
    k=size(Wnew,1);     % 差分向量构建的Hankel矩阵要从k=2维(或di或n维)开始考察
    while 2*k-1<=size(pi_series_dif,2)
        for i=1:k
            pi_delta_FDH(i,:)=pi_series_dif(1,i:k+i-1);
            qi_delta_FDH(i,:)=qi_series_dif(1,i:k+i-1);
        end 
        
     %********** 第3版程序.同时判断p和q两个差分迭代Hankel的秩,取最大者用于计算,可减轻无法算出平均值的状况 *** 
        if (rank(pi_delta_FDH)<k & rank(qi_delta_FDH)<k)                       
            if rank(pi_delta_FDH)>=rank(qi_delta_FDH)
                ri=rank(pi_delta_FDH);             
                delta_i_vec=null(pi_delta_FDH(1:ri+1,1:ri+1),'r')';       % ***注意:写成pi_delta_FDH(1:ri+1,1:ri+1)而不是pi_delta_FDH
                
                % ******注意：有些情况下,这里会算出空向量,这种情况需要处理  (这段代码用于k取n的情况) ********
                    while (isempty(delta_i_vec) && (ri+1)<size(Wnew,1)) || size(delta_i_vec,1)>1                        
                        for rightslideing=1:size(Wnew,1)-ri-1
                            for downsliding=0:rightslideing
                                % 如果主子矩阵计算失败,那么取pi_delta_FDH中右下方的子方阵用于计算
                                delta_i_vec=null(pi_delta_FDH(1+downsliding:ri+1+downsliding,1+rightslideing:ri+1+rightslideing),'r')';
                                if ~isempty(delta_i_vec) && size(delta_i_vec,1)==1
                                    break;
                                end                                  
                            end
                            if ~isempty(delta_i_vec) && size(delta_i_vec,1)==1
                                break;
                            end  
                        end
                        if rightslideing==size(Wnew,1)-ri-1 && downsliding==rightslideing && isempty(delta_i_vec)
                            delta_i_vec=zeros(1,ri+1); 
                        end
                    end
                    if isempty(delta_i_vec) && (ri+1)==size(Wnew,1)
                        delta_i_vec=zeros(1,ri+1);  
                    end
                % *********************************************************************
            else
                ri=rank(qi_delta_FDH);             
                delta_i_vec=null(qi_delta_FDH(1:ri+1,1:ri+1),'r')';    
                
                % ******注意：有些情况下,这里会算出空向量,这种情况需要处理  (这段代码用于k取n的情况) *******
                    while (isempty(delta_i_vec) && (ri+1)<size(Wnew,1)) || size(delta_i_vec,1)>1                        
                        for rightslideing=1:size(Wnew,1)-ri-1
                            for downsliding=0:rightslideing
                                delta_i_vec=null(qi_delta_FDH(1+downsliding:ri+1+downsliding,1+rightslideing:ri+1+rightslideing),'r')';
                                if ~isempty(delta_i_vec) && size(delta_i_vec,1)==1
                                    break;
                                end                                  
                            end
                            if ~isempty(delta_i_vec) && size(delta_i_vec,1)==1
                                break;
                            end  
                        end
                        if rightslideing==size(Wnew,1)-ri-1 && downsliding==rightslideing && isempty(delta_i_vec)
                            delta_i_vec=zeros(1,ri+1); 
                        end
                    end
                    if isempty(delta_i_vec) && (ri+1)==size(Wnew,1)
                        delta_i_vec=zeros(1,ri+1);  
                    end                
                % ***************************************************************************
            end          
            break;         
        else
            pi_delta_FDH=[];      qi_delta_FDH=[];            
            k=k+1; 
        end  
    %***************************************************************************************************   
    end    
    average_value_at_vi=(delta_i_vec*pi_series(1:size(delta_i_vec,2))')/(delta_i_vec*qi_series(1:size(delta_i_vec,2))');  
    
        % ******** 接下来,专门测试在节点vi处计算平均共识值所用的时间 ********
            FDH_pi_temp=[];FDH_qi_temp=[];delta_temp=[]; 
                for u=1:size(Wnew,1)
                    if u==1
                        maxtimestep_vec=0;                    
                    elseif u<ri+1
                        t0=tic;
                        for i=1:u
                            FDH_pi_temp(i,:)=pi_series_dif(1,i:u+i-1);
                            FDH_qi_temp(i,:)=qi_series_dif(1,i:u+i-1);
                        end
                        rankFDHtemp_pi=rank(FDH_pi_temp);
                        rankFDHtemp_qi=rank(FDH_qi_temp);
                        maxtimestep_vec=[maxtimestep_vec,toc(t0)];
                    elseif u==ri+1
                        t0=tic;
                        for i=1:(ri+1)
                            FDH_pi_temp(i,:)=pi_series_dif(1,i:(ri+1)+i-1);
                            FDH_qi_temp(i,:)=qi_series_dif(1,i:(ri+1)+i-1);
                        end
                        rankFDHtemp_pi=rank(FDH_pi_temp);
                        rankFDHtemp_qi=rank(FDH_qi_temp);                        
                        if rankFDHtemp_pi>=rankFDHtemp_qi
                            delta_temp=null(pi_delta_FDH(1:ri+1,1:ri+1),'r')';
                            average_value_temp=(delta_temp*pi_series(1:size(delta_temp,2))')/(delta_temp*qi_series(1:size(delta_temp,2))'); 
                        else
                            delta_temp=null(qi_delta_FDH(1:ri+1,1:ri+1),'r')';
                            average_value_temp=(delta_temp*pi_series(1:size(delta_temp,2))')/(delta_temp*qi_series(1:size(delta_temp,2))');
                        end
                        maxtimestep_vec=[maxtimestep_vec,toc(t0)];                    
                    else
                        maxtimestep_vec=[maxtimestep_vec,0];
                    end               
                    FDH_pi_temp=[];FDH_qi_temp=[];delta_temp=[]; 
                end                
                maxtimestep=max(maxtimestep_vec);
         % *******************************************************************
end