function [di_real,ri,average_value_at_vi,maxtimestep]=computing_Average_consensus(Wnew,eTi,pi_series,qi_series)

%********************** ��vi�������׼��di,delta_i���� **********************
    O_delta_matrix=[];    maxtimestep_vec=[];
    miu_new=sum(Wnew(:,1));                                % ���к�
    t=1;    t_terminal=2*size(Wnew,1)+1;                   % ��"+1"ǰ���õ����ܴ���
    O_delta_matrix(1,:)=eTi*(Wnew-miu_new*eye(size(Wnew,1)));
    while true   
        O_delta_matrix=[O_delta_matrix;O_delta_matrix(t,:)*Wnew];
        t=t+1;             % ��ʾ�Ѿ�����t����������(����),��ʵ��������(ʱ��)��t-1    
        if (t==t_terminal)
            break; 
        end
    end
    AAA_temp=[];
    i=2;
    while i<=size(O_delta_matrix,1)
        AAA_temp=O_delta_matrix(1:i,:);
        if (rank(AAA_temp)<size(AAA_temp,1))
    %         delta_i_real=null(AAA_temp','r')';    % �����ڼ����׼��delta_i 
            di_real=rank(AAA_temp');              % ��׼��di,ֱ����Ϊ�������
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
        pi_series_dif(1,i)=pi_series(1,i+1)-miu_new*pi_series(1,i);     % ��p�����i�й���һ���������
    end   
    for i=1:size(qi_series,2)-1
        qi_series_dif(1,i)=qi_series(1,i+1)-miu_new*qi_series(1,i);     % ��p�����i�й���һ���������
    end  
   
    ri=1;               % ������¼First Hankel�������    
    k=size(Wnew,1);     % �������������Hankel����Ҫ��k=2ά(��di��nά)��ʼ����
    while 2*k-1<=size(pi_series_dif,2)
        for i=1:k
            pi_delta_FDH(i,:)=pi_series_dif(1,i:k+i-1);
            qi_delta_FDH(i,:)=qi_series_dif(1,i:k+i-1);
        end 
        
     %********** ��3�����.ͬʱ�ж�p��q������ֵ���Hankel����,ȡ��������ڼ���,�ɼ����޷����ƽ��ֵ��״�� *** 
        if (rank(pi_delta_FDH)<k & rank(qi_delta_FDH)<k)                       
            if rank(pi_delta_FDH)>=rank(qi_delta_FDH)
                ri=rank(pi_delta_FDH);             
                delta_i_vec=null(pi_delta_FDH(1:ri+1,1:ri+1),'r')';       % ***ע��:д��pi_delta_FDH(1:ri+1,1:ri+1)������pi_delta_FDH
                
                % ******ע�⣺��Щ�����,��������������,���������Ҫ����  (��δ�������kȡn�����) ********
                    while (isempty(delta_i_vec) && (ri+1)<size(Wnew,1)) || size(delta_i_vec,1)>1                        
                        for rightslideing=1:size(Wnew,1)-ri-1
                            for downsliding=0:rightslideing
                                % ������Ӿ������ʧ��,��ôȡpi_delta_FDH�����·����ӷ������ڼ���
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
                
                % ******ע�⣺��Щ�����,��������������,���������Ҫ����  (��δ�������kȡn�����) *******
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
    
        % ******** ������,ר�Ų����ڽڵ�vi������ƽ����ʶֵ���õ�ʱ�� ********
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