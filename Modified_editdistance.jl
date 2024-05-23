
#modified edit distance according to "Recurrence analysis of extreme event-like data"
function M_editdistance(t1, t2,Tao,p0)
    N1 = length(t1)
    N2 = length(t2)

    G = zeros(N1 + 1, N2 + 1)

    G[:, 1] = 0:N1
    G[1, :] = 0:N2

    if N1 > 0 && N2 > 0
        for (i, t1_) in enumerate(t1), (j, t2_) in enumerate(t2)
            try
                cost_del = G[i, j+1] + 1
                cost_add = G[i+1, j] + 1
                cost_shift = G[i, j] + 1/(1+exp(-p0*(abs(t1_ - t2_)-Tao)))
                G[i+1, j+1] = min(cost_del,cost_add,cost_shift)
            catch
                G[i+1, j+1] = min(G[i, j+1] + 1, G[i+1, j] + 1,G[i, j] + 1/(1+exp(-p0*(abs(t1_ - t2_)-Tao))))
                # print("cost function wrong")
            end
        end
        d = G[end,end]
        
    else
        # d = abs(N1 - N2) * 1
        d = 999999; 
    end
    return d
end



function M_distancematrix(t_events, t, w, ws,Tao,p0)

    # calculate parameters for cost function
    # p0 = length(t_events) / (t_events[end] - t_events[1]);
    # p2 = 1.0;
    
    # vector of starting points of the moving window
    time_windows = t[1:ws:end-w];
    
    Δt = mean(diff(t))
    
    # initialize results matrix
    D = 999999*ones(length(time_windows),length(time_windows));

    for (i, w_i) in enumerate(time_windows)
        start = findfirst( >=(w_i), t_events)
        stop = findlast( <(w_i+Δt*w), t_events)
        if isnothing(start) || isnothing(stop)
           #print("start:", start, "\n")
           #print("stop:", stop, "\n")
           #print("no points")
           continue
        end
        if stop < start # window without events
           continue
        end
        i1 = start:stop
        #print("col " , i, " - ", t_events[i1], "\n")
        for (j, w_j) in enumerate(time_windows)
            start = findfirst( >=(w_j), t_events);
            stop = findlast( <(w_j+Δt*w), t_events);
            if isnothing(start) || isnothing(stop)
               continue
            end
            if stop < start # window without events
               continue
            end
            i2 = start:stop
            D[i,j] = M_editdistance(t_events[i1] .- w_i, t_events[i2] .- w_j,Tao,p0)#Note here that the events corresponding to  1-year window are placed in the same time period for easy comparison.
            #if i < 3 & j < 20
            #print("raw " , i, ", ", j, " - ", t_events[i1] .- w_i, ", ", t_events[i2] .- w_j, "\n")
            #print("cost (" , i, ",", j, "): ", D[i,j], "\n")
            #end

        end
    end
    return D
end












