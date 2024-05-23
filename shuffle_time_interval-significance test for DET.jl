
#Note that this code fixes p0=0.05 and then uses the optimal Tau to calculate the significance test.

#-----------------------------------------------------------------------------sig test------------------------------------------------
using Distributed

addprocs(3)

@everywhere using DynamicalSystems
@everywhere using OrdinaryDiffEq
@everywhere using Peaks
@everywhere using Statistics
@everywhere using RecurrenceAnalysis
@everywhere using DelayEmbeddings
@everywhere using DelimitedFiles
@everywhere using CSV
@everywhere using DataFrames
@everywhere using SparseArrays
@everywhere using Distributed
@everywhere using RData
@everywhere using Dates
@everywhere using Random
@everywhere using SharedArrays
@everywhere include("Modified_editdistance.jl")
@everywhere using StatsBase
@everywhere using Distributions


# some parameters
Δt = 1; # sampling time
RR = .05; # preselected recurrence rate
w = 365; ws = 30; # window parameters for edit distance

#sig level
quan = 0.95

obj = load("floodevents_morethan1mm_95th_2002_2016.RData")
fldEvents = obj["floodevents_morethan1mm_95th_2002_2016"];

# @everywhere function exp_Shuffle(t_events)

#     dataDiff = sort(diff(t_events))
#     exp_dist = fit(Exponential, dataDiff)

#     return round.(rand(exp_dist, length(t_events)+1))

# end

@everywhere function exp_Shuffle_Tconst(t_events)

    dataDiff = sort(diff(t_events))
    exp_dist = fit(Exponential, dataDiff)

    TimeInterval =  round.(rand(exp_dist, 5000))

    TimeIntervalCom = cumsum(TimeInterval)

    idx = searchsortedfirst(TimeIntervalCom,5267)-1 #5267 is chosen here because the searchsortedfirst function searches for the first index that is greater than or equal to the specified value.
    return(TimeInterval[1:idx])

end


Best_Tau = CSV.read("MT1mm_fixP0_005_bestTau.csv",header=false,DataFrame);
Best_Tau = Best_Tau.Column1;



###############significance test---------------------------------------------------------
significantGrids = SharedVector{Float64}(fill(NaN, size(fldEvents)[2]))
DETall = SharedVector{Float64}(fill(NaN, size(fldEvents)[2]))
LAMall = SharedVector{Float64}(fill(NaN, size(fldEvents)[2]))
Ratioall =  SharedVector{Float64}(fill(NaN, size(fldEvents)[2]))


task2 = @distributed for i in 1:4#size(fldEvents)[2]
    tr = fldEvents[:,i];
    if isnan(Best_Tau[i])
        continue
    end
    if sum(tr)<3  #Less than 3 games will not be considered
        continue
    end

    det50 = NaN*ones(100);
    
    for j in 1:100
        t = collect(0:Δt:5266); 
        events = tr.==1;
        t_events = t[events];
        getValue = exp_Shuffle_Tconst(t_events)
        #Accumulate to obtain the real position coordinates
        shuffleTevents = cumsum(getValue)
        tnew = t; #Keep the timing the same
        # time_windows = tnew[1:ws:end-w];
        # distance matrix based on modified edit distance
        D = M_distancematrix(shuffleTevents, tnew, w, ws,Best_Tau[i],0.05);
        # RP (Note: heatmap draw map from the top left of the original, i.e. it rotl90() )
        try
        	ε = quantile(vec(D)[vec(D).!=999999],RR); # get threshold based on quantile
			# ε = quantile(vec(D),RR);
            R_ed = D.<ε; 
            # p2 = Plots.heatmap(time_windows .+ w*Δt/2, time_windows .+ w*Δt/2, R_ed);
            R_edPars = sparse(R_ed)
            DET = determinism((R_edPars); lmin=2, theiler = 1)
            LAM = laminarity((R_edPars); lmin=2, theiler = 1)
            if LAM==0
           	 	LAM=1e-6
            end
            det50[j] = DET
        catch ex
            # println(i)
            continue
        end
    end

    #Calculate DET_LAM of original sequence
    t = collect(0:Δt:5266); 
    events = tr.==1;
    t_events = t[events];
    # distance matrix based on modified edit distance
    D = M_distancematrix(t_events, t, w, ws,Best_Tau[i],0.05);
    # RP (Note: heatmap draw map from the top left of the original, i.e. it rotl90() )
    try
        ε = quantile(vec(D)[vec(D).!=999999],RR); # get threshold based on quantile
        # ε = quantile(vec(D),RR); # get threshold based on quantile
        R_ed = D.<ε;
        # p2 = Plots.heatmap(time_windows .+ w*Δt/2, time_windows .+ w*Δt/2, R_ed);
        R_edPars = sparse(R_ed)
        DET = determinism((R_edPars); lmin=2, theiler = 1)
        LAM = laminarity((R_edPars); lmin=2, theiler = 1)
        if LAM==0
           	LAM=1e-6
        end
        original_DET =  DET
        DETall[i] = DET
        LAMall[i] = LAM
        Ratioall[i] = DET/LAM
        #
        if size(filter(!isnan, det50))[1]==0
           continue
        else
           qn = quantile(filter(!isnan, det50), quan)
           if original_DET> qn
                significantGrids[i] = DETall[i]
           end

        end
    catch ex
        println("There is a problem with the $i -th grid point!")
        # println(i)
        continue
    end
    
end

# wait for the task finished
wait(task2)

#remember to change the data path
CSV.write("MT1mm_exp_Shuffle_p0_005_bestTau_DET_significantGrids95th-Tconst.csv",Tables.table(significantGrids),writeheader=false)
# CSV.write("/work/home/bjsfdxzjx/JianxinFile/P0_005_findBestTau/result2/exp_Shuffle_p0_005_bestTau_DET.csv",Tables.table(DETall),writeheader=false)
# CSV.write("/work/home/bjsfdxzjx/JianxinFile/P0_005_findBestTau/result2/exp_Shuffle_p0_005_bestTau_LAM.csv",Tables.table(LAMall),writeheader=false)
# CSV.write("/work/home/bjsfdxzjx/JianxinFile/P0_005_findBestTau/result2/exp_Shuffle_p0_005_bestTau_DETratioLAM.csv",Tables.table(Ratioall),writeheader=false)



























