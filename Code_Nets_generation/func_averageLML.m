function [average_loopslength,MaximumLoopLength] = func_averageLML(G,A,p_graph)
%FUNC_AVERAGELML 此处显示有关此函数的摘要
% p_graph=plot(G,'Layout','force','EdgeColor',[0 0.4 0.8],'LineWidth',0.5,'NodeColor','k','NodeFontSize',8);

loops_sum=[];
timefail=1;
% PP = parpool(1);

for i=size(A,2):-1:1

        PF = parfeval(@allcycles,2,G,'MaxNumCycles',1,'MinCycleLength',i);       % 用后台并行运算的方法，找到1个最长的环
        t0=tic;                                                                                 % 设置单次计算失败的时间限制
        while true           
            if toc(t0)>timefail && (strcmp(PF.State,'running') )
                cancel(PF);
                longestloop_temp=[];
                edgelongestloop=[];
                break;                                              % 超出时间限制，就跳出
            elseif strcmp(PF.State,'finished')
                [longestloop_temp,edgelongestloop]=fetchOutputs(PF);
                break;                                              % 不超出时间限制且成功找到最长环，跳出            
            else
                longestloop_temp=[];
                edgelongestloop=[];
            end             
        end
%         [longestloop_temp,edgelongestloop]=allcycles(G,'MaxNumCycles',1,'MinCycleLength',i,'MaxCycleLength',i); % 找到1个最长的环
        LongestLoop_temp_index=cell2mat(longestloop_temp);
        MaximumLoopLength=length(LongestLoop_temp_index);               % 最长环的长度
        if ~isempty(longestloop_temp)
            loops_sum=MaximumLoopLength^2;
            break;
        end
        clear PF t0;
%         delete(PP);
 end
highlight(p_graph,'Edges',edgelongestloop{1},'EdgeColor','r','LineWidth',1.8,'NodeColor','b','MarkerSize',5)   % 红色标出其中一个最长的环


network_index=(1:size(A,2));
LongestLoopMissed_index = find(ismember(network_index,LongestLoop_temp_index)==0);  % 找到哪些节点不在最长环上

for i=LongestLoopMissed_index
    
    A_temp=circshift(A,[-(i-1),-(i-1)]);            % 对原矩阵向下和向右平移(i-1)位
    G_temp = digraph(A_temp');

    for j=MaximumLoopLength:-1:1
        PF = parfeval(@allcycles,1,G_temp,'MaxNumCycles',1,'MinCycleLength',j,'MaxCycleLength',j);
%         PF = parfeval(backgroundPool,@allcycles,1,G_temp,'MaxNumCycles',1,'MinCycleLength',j);   % 有时会找到更长的环,这样的平均环长可能大于最长环长
        t0=tic;
        while true           
            if toc(t0)>timefail && strcmp(PF.State,'running')
                cancel(PF);
                cycle_temp=[];
                break;
            elseif strcmp(PF.State,'finished')
                cycle_temp=fetchOutputs(PF);
                break;
            else
                cycle_temp=[];
            end             
        end
        Cycle_temp_index=cell2mat(cycle_temp);
        if ~isempty(cycle_temp) && Cycle_temp_index(1)==1
            loops_sum=loops_sum+length(Cycle_temp_index);
            break;
        end
        clear PF t0;
    end

end
average_loopslength=loops_sum/size(A,2);

end

