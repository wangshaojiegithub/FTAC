function [di_real,ri_max,ri_max_source_vec,LeftVecUni_real_at_vi,LeftVecUni_at_vi,maxtimestep]=computing_LeftEigenvector_from_channels(vi,eTi,W)

    t_terminal=2*size(W,1)+1;                   % 在"+1"前设置迭代总次数
    p_matrix=[];  q_matrix=[]; O_alpha_matrix=[]; O_delta_matrix=[];    maxtimestep_vec=[];

    miu=sum(W(1,:));                    % 求行和
    I=eye(size(W,1));
    for vi=1:size(W,1)    
        eval(['p',num2str(vi),'_vec=I(:,vi)']);
    end 
    q0_vec([1:size(W,1)],1)=1;

    %********** 开始计算：在vi处计算标准的di, delta_i向量和left vec(表示成left_vi_real_Uni)**********
    t=1;
    O_alpha_matrix(1,:)=eTi;
    O_delta_matrix(1,:)=eTi*(W-miu*eye(size(W,1)));
    while true
        O_alpha_matrix=[O_alpha_matrix;O_alpha_matrix(t,:)*W];
        O_delta_matrix=[O_delta_matrix;O_delta_matrix(t,:)*W];
        t=t+1;             % 表示已经具有t个连续数据(向量),其实迭代次数(时间)是t-1    
        if (t==t_terminal)
            break; 
        end
    end 
    
    % 计算标准的delta_i向量(与初值无关)
    AAA_temp=[];
    i=2;
    while true
        while i<=size(O_delta_matrix,1)
            AAA_temp=O_delta_matrix(1:i,:);
            if (rank(AAA_temp)<size(AAA_temp,1))
    %             delta_i_real=null(AAA_temp','r')';    % 可用于计算标准的delta_i       
    
                % ******* 换成数值稳定性更好的null函数 ******************
                delta_i_real=null(AAA_temp')';  
                if ~isempty(delta_i_real)
                NormalizingNumber=1/delta_i_real(end);  
                delta_i_real=delta_i_real*NormalizingNumber;
                end
                % ***************************************************
                
                di_real=rank(AAA_temp');              % 标准的di,直接作为输出数据
                break;
            else
                AAA_temp=[];
                i=i+1;
            end        
        end
        % *****在vi处计算出的标准的左特征向量(与初值无关),可作为接下来用观测序列计算出的左特征向量的对照值 *****
        leftvec_real_at_vi=delta_i_real*O_alpha_matrix(1:size(delta_i_real,2),:); 
    
        if size(leftvec_real_at_vi,1)==size(W,1)
            LeftVecUni_real_at_vi=leftvec_real_at_vi/sum(leftvec_real_at_vi);  
            break;
        else
            leftvec_real_at_vi=[];
            AAA_temp=[];
            i=i+1;
        end
    end
    
        %*************** 求标准的左特征向量,可替换为如下代码 *************************
%             eigenvalues=eig(W);eigenvalues=roundn(eigenvalues,-4);
%             indexeigen=find(eigenvalues==miu);
%             [eigenvec,~]=eig(W');
%             LeftVecUni_real_at_vi=eigenvec(:,indexeigen)/sum(eigenvec(:,indexeigen));
        %*********************************************
        
    % ****************************************************************************************

    % 在vi处,接下来的每个ci表示第i个channel(共n个)
    for  cj=1:size(W,1)    
        eval(['pi_vec=p',num2str(cj),'_vec']);
        t=1;   t_terminal=2*size(W,1)+1;
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
        pi_series=eTi*p_matrix;
        qi_series=eTi*q_matrix;

        pi_series_dif=[];
        for i=1:size(pi_series,2)-1
            pi_series_dif(1,i)=pi_series(1,i+1)-miu*pi_series(1,i);% 由p矩阵第i行构建一个差分向量
        end    
        pi_delta_FDH=[];    % 差分Hankel矩阵
        ri_j=1;             % 用来记录vi处计算第j通道值的 差分FDH矩阵的秩

        k=size(W,1);            % *********差分向量构建的Hankel矩阵要从k=2维(或di或n维)开始考察
        while 2*k-1<=size(pi_series_dif,2)
            for i=1:k
                pi_delta_FDH(i,:)=pi_series_dif(1,i:k+i-1);
            end 
            if (rank(pi_delta_FDH)<k)
                break; 
            else
                pi_delta_FDH=[];    % 注意如果该Hankel不是缺陷的，要在下一次循环清空，以免矩阵维数被固定
                k=k+1;
            end   
        end       
        ri_j=rank(pi_delta_FDH);   % 其实现在的Hankel矩阵的秩是r=k-1。
%         delta_i_j=null(pi_delta_FDH(1:ri_j+1,1:ri_j+1),'r')';      
        
        % *************  数值稳定性更好的 null函数 ******************
        delta_i_j=null(pi_delta_FDH(1:ri_j+1,1:ri_j+1))';
        if ~isempty(delta_i_j)
        NormalizingNumber=1/delta_i_j(end);  
        delta_i_j=delta_i_j*NormalizingNumber;
        end
        % *********************************************************
        
        % ****** 注意：有些情况下,这里会算出空向量,这种情况需要处理 (这段代码用于k取n的情况)*****************
            while (isempty(delta_i_j) && (ri_j+1)<size(W,1)) || size(delta_i_j,1)>1                        
                for rightslideing=1:size(W,1)-ri_j-1
                    for downsliding=0:rightslideing
                    % 如果主子矩阵计算失败,那么取pi_delta_FDH中右下方的子方阵用于计算
%                         delta_i_j=null(pi_delta_FDH(1+downsliding:ri_j+1+downsliding,1+rightslideing:ri_j+1+rightslideing),'r')';

                        % *************  数值稳定性更好的 null函数 ******************
                        delta_i_j=null(pi_delta_FDH(1+downsliding:ri_j+1+downsliding,1+rightslideing:ri_j+1+rightslideing))';
                        if ~isempty(delta_i_j)
                        NormalizingNumber=1/delta_i_j(end);  
                        delta_i_j=delta_i_j*NormalizingNumber;
                        end
                        % *********************************************************
                        
                        if ~isempty(delta_i_j) && size(delta_i_j,1)==1
                            break;
                        end                                  
                    end
                    if ~isempty(delta_i_j) && size(delta_i_j,1)==1
                        break;
                    end  
                end
                if rightslideing==size(W,1)-ri_j-1 && downsliding==rightslideing && isempty(delta_i_j)
                    delta_i_j=zeros(1,ri_j+1); 
                end
            end
            if isempty(delta_i_j) && (ri_j+1)==size(W,1)
                delta_i_j=zeros(1,ri_j+1);  
            end         
         % *****************************************************************
         
        left_i_j=(delta_i_j*pi_series(1:size(delta_i_j,2))')/(delta_i_j*qi_series(1:size(delta_i_j,2))'); % 标准的(归一化)左特征向量第i个分量

        % ******************************* 接下来,专门测试在节点vi处计算左特征向量的分量xj所用的时间 *******************
                FDH_temp=[];delta_temp=[];      
                if cj==1
                    for u=1:size(W,1)
                        if u==1
                            maxtimestep_vec=0;                    
                        elseif u<ri_j+1
                            t0=tic;
                            for i=1:u
                                FDH_temp(i,:)=pi_series_dif(1,i:u+i-1);
                            end
                            rankFDHtemp=rank(FDH_temp);
                            maxtimestep_vec=[maxtimestep_vec,toc(t0)];
                        elseif u==ri_j+1
                            t0=tic;
                            for i=1:(ri_j+1)
                                FDH_temp(i,:)=pi_series_dif(1,i:(ri_j+1)+i-1);
                            end
                            rankFDHtemp=rank(FDH_temp);
%                             delta_temp=null(pi_delta_FDH(1:ri_j+1,1:ri_j+1),'r')';
                            
                            % *************  数值稳定性更好的 null函数 ******************
                            delta_temp=null(pi_delta_FDH(1:ri_j+1,1:ri_j+1))';
                            if ~isempty(delta_temp)
                            NormalizingNumber=1/delta_temp(end);  
                            delta_temp=delta_temp*NormalizingNumber;
                            end
                            % *********************************************************
                        
                            left_temp=(delta_temp*pi_series(1:size(delta_temp,2))')/(delta_temp*qi_series(1:size(delta_temp,2))');
                            maxtimestep_vec=[maxtimestep_vec,toc(t0)];                    
                        else
                            maxtimestep_vec=[maxtimestep_vec,0];
                        end               
                        FDH_temp=[];delta_temp=[];
                    end                
                 else
                    for u=1:size(W,1)
                        if u==1
                            maxtimestep_vec(u)=0;                    
                        elseif u<ri_j+1
                            t0=tic;
                            for i=1:u
                                FDH_temp(i,:)=pi_series_dif(1,i:u+i-1);
                            end
                            rankFDHtemp=rank(FDH_temp);
                            maxtimestep_vec(u)=maxtimestep_vec(u)+toc(t0);
                        elseif u==ri_j+1
                            t0=tic;
                            for i=1:(ri_j+1)
                                FDH_temp(i,:)=pi_series_dif(1,i:(ri_j+1)+i-1);
                            end
                            rankFDHtemp=rank(FDH_temp);
%                             delta_temp=null(pi_delta_FDH(1:ri_j+1,1:ri_j+1),'r')';
                            
                            % *************  数值稳定性更好的 null函数 ******************
                            delta_temp=null(pi_delta_FDH(1:ri_j+1,1:ri_j+1))';
                            if ~isempty(delta_temp)
                            NormalizingNumber=1/delta_temp(end);  
                            delta_temp=delta_temp*NormalizingNumber;
                            end
                            % *********************************************************
                            
                            left_temp=(delta_temp*pi_series(1:size(delta_temp,2))')/(delta_temp*qi_series(1:size(delta_temp,2))');
                            maxtimestep_vec(u)=maxtimestep_vec(u)+toc(t0);                  
                        else
                            maxtimestep_vec(u)=maxtimestep_vec(u)+0;
                        end   
                        FDH_temp=[];delta_temp=[];
                    end                
                end         
        % ***************************************************************************************
        %   记录要输出的数据
        if cj==1
            ri_max=ri_j;
            ri_max_source_vec=ri_j;        
            LeftVecUni_at_vi=left_i_j;
        else
            ri_max=max([ri_max,ri_j]);
            ri_max_source_vec=[ri_max_source_vec,ri_j];        
            LeftVecUni_at_vi=[LeftVecUni_at_vi,left_i_j];
        end      
        p_matrix=[];  q_matrix=[]; O_alpha_matrix=[]; O_delta_matrix=[];     delta_i_j=[];
    end 
        maxtimestep=max(maxtimestep_vec);

end

