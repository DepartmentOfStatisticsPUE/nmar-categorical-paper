## riddles et al nmar

using Distributions
using Random
using DataFrames
using StatsBase
using Statistics
using FreqTables
using DataFramesMeta

## population data
Random.seed!(123);
strata_names = [11,12,21,22]
strata_sizes = 500 .* [50,100,150,200]
## x = vcat(fill.(strata_names, strata_sizes)...)
x = inverse_rle(strata_names, strata_sizes)
df = DataFrame(x = x)
df.x1 = ifelse.(SubString.(string.(df.x), 1, 1) .== "1", 1, 2)
df.x2 = ifelse.(SubString.(string.(df.x), 2, 2) .== "1", 1, 2)
df.y = vcat(
    wsample([1,2], [0.7,0.3], strata_sizes[1]),  
    wsample([1,2], [0.5,0.5], strata_sizes[2]), 
    wsample([1,2], [0.3,0.7], strata_sizes[3]),  
    wsample([1,2], [0.4,0.6], strata_sizes[4])
)
df.eta = -0.4 .* (df.x1 .== 2) .- 0.8 .* (df.y .== 2)
df.rho = 1 ./ (1 .+ exp.(df.eta))

## sampling Random.seed!(b)

B = 500
sim_res = zeros(B,2)
for b in 1:B
    Random.seed!(b);
    df.resp = [rand(Bernoulli(i)) for i in df.rho]
    df_sampl = by(df, [:resp, :x1, :x2, :y], n = :resp => length)
    
    df_sampl_obs = df_sampl[df_sampl.resp .== 1, :]
    df_sampl_obs = @transform(groupby(df_sampl_obs, [:x1, :x2]), p_hat = :n/sum(:n))
    
    df_sampl_obs.p_hat = by(df_sampl_obs, [:x1, :x2], p_hat = :n => (x -> x/sum(x))).p_hat
    df_sampl_nonobs = by(df_sampl[df_sampl.resp .== 0,:], [:x1, :x2], m = :n => sum)
    
    df_sampl_obs = leftjoin(df_sampl_obs, df_sampl_nonobs, on = [:x1, :x2])
    df_sampl_obs.O = 1
    
    O_start = df_sampl_obs.O
    
    ## interations
    
    for iter in 1:10000
        df_sampl_obs = @transform(groupby(df_sampl_obs, [:x1, :x2]), m_hat = :m .* :p_hat .* :O / sum(:p_hat .* :O))
        df_sampl_obs = @transform(groupby(df_sampl_obs, [:y, :x1]), O = sum(:m_hat) / sum(:n))
        dif = sum((O_start - df_sampl_obs.O).^2)
        
        if (dif < sqrt(eps()))
            println("Converged on interation: ", iter)
            break
        end
        O_start = df_sampl_obs.O
    end
    
    result = by(sort!(df_sampl_obs,:y), :y, x -> sum(x.n .* (1 .+ x.O))) 
    sim_res[b, 1] =  result.x1[1]/(result.x1[1]+result.x1[2])
    sim_res[b, 2] =  result.x1[2]/(result.x1[1]+result.x1[2])
end 

npar_bias = mean(sim_res, dims = 1) - prop(convert(Array, freqtable(df.y))')
npar_var = var(sim_res, dims = 1)
npar_mse = npar_bias.^2 + npar_var