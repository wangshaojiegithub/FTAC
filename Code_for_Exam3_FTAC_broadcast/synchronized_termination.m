function [S1_vi_ready_at_t_vec,S1_vi_stop_at_t_vec]=synchronized_termination(W,dist,epsi_max,S1_ri_max_vec)

t=0;          t_final=40;                   % 注意：t是从0开始
count_vi_vec_pre  =    zeros(size(W,1),1);     % 上一时刻的count_vi_vec(设置它的意义是为了避免在vi循环中错误地利用当前时刻已更新的邻居的count值)
count_vi_vec   =       zeros(size(W,1),1);     % 当前时刻的count_vi_vec
S1_vi_ready_at_t_vec = zeros(size(W,1),1);
S1_vi_stop_at_t_vec =  zeros(size(W,1),1);

while (t~=t_final)
    for  vi=1:size(W,1)
        if t<2*(S1_ri_max_vec(vi)+1)-1 & t>=epsi_max(vi)            
            count_vi_vec(vi)=epsi_max(vi);     % 给count_vi赋初值
        elseif t>=2*(S1_ri_max_vec(vi)+1)-1 & count_vi_vec(vi)~=-1
            if t==2*(S1_ri_max_vec(vi)+1)-1
                S1_vi_ready_at_t_vec(vi)=t;        % 记录vi准备停止的时间,就是完成本地计算的时刻t=2(ri_max+1)-1
            end 
            neighbors_index=find(dist(vi,:)==1 | dist(vi,:)==0);        % vi及其内邻居的索引(是一个行向量)
            count_vi_vec(vi)=max(count_vi_vec_pre(neighbors_index))-1;  % 给count_vi赋新值            
            if (count_vi_vec(vi)==-1 & S1_vi_stop_at_t_vec(vi)==0)
                S1_vi_stop_at_t_vec(vi)=t;                  % 记录vi停止的时刻
            end
        end  
    end   
    count_vi_vec_pre=count_vi_vec;  % 将当前时刻的count值存储在count_vi_vec_pre中
    t=t+1;    
end




end
