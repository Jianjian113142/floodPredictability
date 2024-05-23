using Distributed

addprocs(64)


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

# some parameters
Δt = 1; # sampling time
RR = .05; # preselected recurrence rate
w = 365; ws = 30; # window parameters for edit distance


obj = load("floodevents_morethan1mm_95th_2002_2016.RData")
fldEvents = obj["floodevents_morethan1mm_95th_2002_2016"];

RatioAll = SharedVector{Float64}(fill(NaN, size(fldEvents)[2]))
Tao_bset = SharedVector{Float64}(fill(NaN, size(fldEvents)[2]))
DETall = SharedVector{Float64}(fill(NaN, size(fldEvents)[2]))
LAMall = SharedVector{Float64}(fill(NaN, size(fldEvents)[2]))


TAO = CSV.read("globalfld_Tau_meanTimeInterval.csv",header=false,DataFrame);
TAO = TAO.Column1;



task1 = @distributed for i in 1:size(fldEvents)[2]
    # # time axis 
    tr = fldEvents[:,i];
    t = collect(0:Δt:5266); 

    # create event series
    events = tr.==1; # indices and values of local maxima
    t_events = t[events]; # time points of local maxima

    Temp = NaN*ones(size(TAO)[1]);
	DETtemp = NaN*ones(size(TAO)[1]);
    LAMtemp = NaN*ones(size(TAO)[1]);
    #set Tau range
    if isnan(TAO[i])
        continue
    end
    Taurange = 0:1:TAO[i]
   
    for n in 1:length(Taurange)
            # distance matrix based on modified edit distance
            D = M_distancematrix(t_events, t, w, ws,Taurange[n],0.05);
            # RP (Note: heatmap draw map from the top left of the original, i.e. it rotl90() )
            try
                ε = quantile(vec(D)[vec(D).!=999999],RR); # get threshold based on quantile
                R_ed = D.<ε;
                # p2 = Plots.heatmap(time_windows .+ w*Δt/2, time_windows .+ w*Δt/2, R_ed);
                R_edPars = sparse(R_ed)
                DET = determinism((R_edPars); lmin=2, theiler = 1)
                LAM = laminarity((R_edPars); lmin=2, theiler = 1)
                if LAM==0
       					 LAM = 1e-6
   					  end
                Temp[n] = DET/LAM
                DETtemp[n] = DET
                LAMtemp[n] = LAM
            catch ex
                println("There is no events in the $i -th grid point!")
                # println(i)
                continue
            end
        
    end
    

    if size(filter(!isnan, Temp))[1]==0
        RatioAll[i] = NaN
        Tao_bset[i] = NaN
        #p0_bset[i] = NaN
        DETall[i] = NaN
        LAMall[i] = NaN
    elseif size(filter(isnan, Temp))[1]>0
        Temp[isnan.(Temp)] .= -1
        RatioAll[i] = maximum(Temp);
        Tao_bset[i] = Taurange[argmax(Temp)]
        DETall[i] = DETtemp[argmax(Temp)]
        LAMall[i] = LAMtemp[argmax(Temp)]
    else
        RatioAll[i] = maximum(Temp);
        Tao_bset[i] = Taurange[argmax(Temp)]
        DETall[i] = DETtemp[argmax(Temp)]
        LAMall[i] = LAMtemp[argmax(Temp)]
	end
end

# wait for the task finished
wait(task1)

CSV.write("MT1mm_fixP0_005_DET_ratio_LAM.csv",Tables.table(RatioAll),writeheader=false)
CSV.write("MT1mm_fixP0_005_bestTau.csv",Tables.table(Tao_bset),writeheader=false)
CSV.write("MT1mm_fixP0_005_DET.csv",Tables.table(DETall),writeheader=false)
CSV.write("MT1mm_fixP0_005_LAM.csv",Tables.table(LAMall),writeheader=false)





