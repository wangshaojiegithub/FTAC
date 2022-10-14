function [S1_vi_ready_asyn_vec,S1_vi_stop_asyn_vec]=asynchronized_termination(W,dist,S1_di_vec)

t=0;          t_final=40;                   % 注意：t是从0开始

S1_vi_ready_asyn_vec = zeros(size(W,1),1);
S1_vi_stop_asyn_vec =  zeros(size(W,1),1);

while (t~=t_final)
    for  vi=1:size(W,1)
        if t==2*(S1_di_vec(vi)+1)-1            
            S1_vi_ready_asyn_vec(vi)=t;        % 记录vi准备停止的时间,就是完成本地计算的时刻t=2(di+1)-1
        elseif t>2*(S1_di_vec(vi)+1)-1 & t-S1_vi_ready_asyn_vec(vi)==2*(S1_di_vec(vi)+1)-1
            S1_vi_stop_asyn_vec(vi)=t;                  % 记录vi停止的时刻
            
        end  
    end   
    t=t+1;    
end


end
