using DataFrames, CSV, StatsBase

data_out = CSV.File("/home/arosen/Documents/hurdleModelExplore/data/ScaleWithOutcomesOSF.csv") |> DataFrame
names(data_out)

response = data_out[:, 2:15] # Select 14 item responses

tabs_all = map(eachcol(response) .=> table)
tabs_all = map(tabs_all) .|> proptable
tabs_all

patterns = combine(groupby(response, [:x1, :x2, :x3, :x4, :x5, :x6, :x7, :x8, :x9, :x10, :x11, :x12, :x13, :x14]), n = nrow)
patterns = sort!(patterns, [:x1, :x2, :x3, :x4, :x5, :x6, :x7, :x8, :x9, :x10, :x11, :x12, :x13, :x14])
rename!(patterns, :n => :r)

x = Matrix(patterns[:, 1:14])
r = patterns.r

nitems = size(x, 2)
n = sum(r)
data = DataFrame(I(x), r = r) # Dataframe with all item responses and frequencies

theta1 = -6:0.25:6 # Quadrature points for theta1
theta2 = -6:0.25:6 # Quadrature points for theta2
theta = collect(Iterators.product(theta1, theta2))
theta = DataFrame(theta, :auto) |> rename!(_, [:theta1, :theta2])

nitemsgrm = nitems
ncatgrm = maximum(x) + 1
ncat2PL = 2
mat = zeros(nitems, length(theta.theta1))
itemtraceGRM = [copy(mat) for _ in 1:ncatgrm]
itemtrace2PL = [copy(mat) for _ in 1:ncat2PL]
itemtrace = [copy(mat) for _ in 1:(ncatgrm + 1)]

function trace_line_pts_2PL(a_z, b_z, theta)
    for j in 1:nitems
        for k in 0:(ncatgrm - 1)
            if k == 0
                itemtrace2PL[1][j, :] .= 1 .- exp.(a_z[j, :] .* (theta[:, 1] .- b_z[j, 1])) ./ (1 .+ exp.(a_z[j, :] .* (theta[:, 1] .- b_z[j, 1])))
            else
                itemtrace2PL[2][j, :] .= exp.(a_z[j, :] .* (theta[:, 1] .- b_z[j, 1])) ./ (1 .+ exp.(a_z[j, :] .* (theta[:, 1] .- b_z[j, 1])))
            end
        end
    end
    return itemtrace2PL
end

function trace_line_pts_grm(a, b, theta)
    for j in 1:nitems
        for k in 0:(ncatgrm - 1)
            if k == 0
                itemtraceGRM[k + 1][j, :] .= 1 .- exp.(a[j, :] .* (theta[:, 2] .- b[j, k + 1])) ./ (1 .+ exp.(a[j, :] .* (theta[:, 2] .- b[j, k + 1])))
            elseif k == 1
                itemtraceGRM[k + 1][j, :] .= exp.(a[j, :] .* (theta[:, 2] .- b[j, k])) ./ (1 .+ exp.(a[j, :] .* (theta[:, 2] .- b[j, k]))) .- exp.(a[j, :] .* (theta[:, 2] .- b[j, k + 1])) ./ (1 .+ exp.(a[j, :] .* (theta[:, 2] .- b[j, k + 1])))
            elseif k == 2
                itemtraceGRM[k + 1][j, :] .= exp.(a[j, :] .* (theta[:, 2] .- b[j, k])) ./ (1 .+ exp.(a[j, :] .* (theta[:, 2] .- b[j, k]))) .- exp.(a[j, :] .* (theta[:, 2] .- b[j, k + 1])) ./ (1 .+ exp.(a[j, :] .* (theta[:, 2] .- b[j, k + 1])))
            elseif k == 3
                itemtraceGRM[k + 1][j, :] .= exp.(a[j, :] .* (theta[:, 2] .- b[j, k])) ./ (1 .+ exp.(a[j, :] .* (theta[:, 2] .- b[j, k])))
            end
        end
    end
    return itemtraceGRM
end

function trace_line_pts(a, b, a_z, b_z, theta)
    for j in 1:nitems
        for k in 0:(ncatgrm - 1)
            if k == 0
                itemtrace[k + 1][j, :] .= (1 .- trace_line_pts_2PL(a_z, b_z, theta)[2][j, :]) .+ 
                    (trace_line_pts_2PL(a_z, b_z, theta)[2][j, :]) .* trace_line_pts_grm(a, b, theta)[k + 1][j, :]
            elseif k == 1 || k == 2 || k == 3
                itemtrace[k + 1][j, :] .= (trace_line_pts_2PL(a_z, b_z, theta)[2][j, :]) .* trace_line_pts_grm(a, b, theta)[k + 1][j, :]
            end
        end
    end
    return itemtrace
end

