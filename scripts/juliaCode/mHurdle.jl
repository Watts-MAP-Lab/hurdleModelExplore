using Combinatorics, Distributions, CSV, DataFrames, Optim, LinearAlgebra

data_out = CSV.File("./data/ScaleWithOutcomesOSF.csv") |> DataFrame


function trace_line_pts_grm(a, b, theta)
    itemtraceGRM = zeros(size(a)[1], size(theta)[1], 3)
    for j in 1:size(a, 1)
        for k in 1:3
            if k == 1
                itemtraceGRM[j, :,1] .= 1 .- exp.(a[j, :] .* (theta[:, 2] .- b[j, k])) ./ (1 .+ exp.(a[j, :] .* (theta[:, 2] .- b[j, k])))
            elseif k == 2
                itemtraceGRM[j, :,2] .= exp.(a[j, :] .* (theta[:, 2] .- b[j, k-1])) ./ (1 .+ exp.(a[j, :] .* (theta[:, 2] .- b[j, k-1]))) .- exp.(a[j, :] .* (theta[:, 2] .- b[j, k])) ./ (1 .+ exp.(a[j, :] .* (theta[:, 2] .- b[j, k])))
            elseif k == 3
                itemtraceGRM[j, :,3] .= exp.(a[j, :] .* (theta[:, 2] .- b[j, k-1])) ./ (1 .+ exp.(a[j, :] .* (theta[:, 2] .- b[j, k-1])))
            end
        end
    end
    return itemtraceGRM
end

function trace_line_pts_2PL(a_z, b_z, theta)
    nitems = size(a_z, 1);
    #itemtrace2PL = [zeros(size(a_z, 1), size(theta, 1)), zeros(size(a_z, 1), size(theta, 1))]
    itemtrace2PL = zeros(nitems, size(theta, 1), 2);
    for j in 1:nitems
        for k in 0:3
            if k == 0
                itemtrace2PL[j, :,1] .= 1 .- exp.(a_z[j, :] .* (theta[:, 1] .- b_z[j, 1])) ./ (1 .+ exp.(a_z[j, :] .* (theta[:, 1] .- b_z[j, 1])))
            else
                itemtrace2PL[j, :,2] .= exp.(a_z[j, :] .* (theta[:, 1] .- b_z[j, 1])) ./ (1 .+ exp.(a_z[j, :] .* (theta[:, 1] .- b_z[j, 1])))
            end
        end
    end
    return itemtrace2PL
end


function trace_line_pts(a, b, a_z, b_z, theta)
    itemtrace = zeros(14,size(theta)[1],4)
    for j in 1:size(a, 1)
        for k in 0:3
            if k == 0
                itemtrace[j, :,k+1] .= 1 .- trace_line_pts_2PL(a_z, b_z, theta)[j, :,2]
            elseif k == 1
                itemtrace[j, :,k+1] .= trace_line_pts_2PL(a_z, b_z, theta)[j, :,2] .* trace_line_pts_grm(a, b, theta)[j, :,1]
            elseif k == 2
                itemtrace[j, :,k+1] .= trace_line_pts_2PL(a_z, b_z, theta)[j, :,2] .* trace_line_pts_grm(a, b, theta)[j, :,2]
            elseif k == 3
                itemtrace[j, :,k+1] .= trace_line_pts_2PL(a_z, b_z, theta)[j, :,2] .* trace_line_pts_grm(a, b, theta)[j, :,3]
            end
        end
    end
    return itemtrace
end

function multivariate_normal_pdf(x, μ, Σ)
    """
    Calculates the probability density function (PDF) of a multivariate normal distribution.

    Args:
        x: A vector representing a data point.
        μ: A vector representing the mean of the distribution.
        Σ: A matrix representing the covariance matrix of the distribution.

    Returns:
        The PDF value at the given point.
    """

    n = size(x)[1]
    if n != length(μ) || n != size(Σ)[1] || n != size(Σ)[2]
        error("Input dimensions must match")
    end

    # Ensure Σ is positive definite
    if !isposdef(Σ)
        error("Covariance matrix must be positive definite")
    end

    # Calculate the determinant and inverse of Σ
    det_Σ = det(Σ)
    inv_Σ = inv(Σ)

    # Compute the exponent term
    exponent = -0.5 * (x .- μ)' * inv_Σ * (x - μ)

    # Calculate the PDF
    pdf = (1 / sqrt(det_Σ * (2π)^n)) * exp(exponent)

    return pdf
end

## Now try to get these to run
a = ones(14)

b = [[1,1,1,1,1,1,1,1,1,1,1,1,1,1] [2,2,2,2,2,2,2,2,2,2,2,2,2,2]]

a_z = ones(14)

b_z = ones(14).*-1

theta1 = -6:0.25:6  # Quadrature points for theta1
theta2 = -6:0.25:6
theta = collect(Combinatorics.permutations(theta1,2))
d = collect(Iterators.product(theta1, theta2))
d = reshape(d, 2401)
theta = [d[i][j] for i in 1:length(d), j in 1:length(d[1])]
trace_line_pts(a, b, a_z, b_z, theta)

a = ones(14)

b = [[1,1,1,1,1,1,1,1,1,1,1,1,1,1] [2,2,2,2,2,2,2,2,2,2,2,2,2,2]]

a_z = ones(14)

b_z = ones(14).*-1
rho_exp = 0.4
nitems = size(testitems)[2]
nParmsPerItemGRM = 3
nParmsPerItem2PL = 2
p = zeros(nitems*nParmsPerItemGRM + nitems*nParmsPerItem2PL + 1); 
nitemsgrm = 14
ncatgrm = 3
for j in 1:nitems 
  p[(j-1)*nParmsPerItemGRM + 1] = a[j]
  p[nitemsgrm*nParmsPerItemGRM + (j-1)*nParmsPerItem2PL + 1] = a_z[j]
  p[nitemsgrm*nParmsPerItemGRM + (j-1)*nParmsPerItem2PL + 2] = b_z[j]
  for k in 1:(ncatgrm-1)
    p[(j-1)*nParmsPerItemGRM + 1 + k] = b[j,k]
  end 
end
p[nitems*nParmsPerItemGRM + nitems*nParmsPerItem2PL + 1] = rho_exp

## Now grab the cross tabs for all response patterns

## Now figure out the LL function here
function ll_grm_ip(p, testdata, theta)
    nParmsPerItemGRM = 3
    nParmsPerItem2PL = 2
    ncatgrm = 3
    r = combine(groupby(testdata, [:x1, :x2, :x3, :x4, :x5, :x6, :x7, :x8, :x9, :x10, :x11, :x12, :x13, :x14]), nrow)
    nitemsgrm=nitems = size(testdata)[2]
    a = fill(-1.0, nitems, 1)
    b = fill(-1.0, nitems, 4 - 1)
    a_z = fill(-1.0, nitems, 1)
    b_z = fill(-1.0, nitems, 1)
    rho_exp = 0.5

    ## Fix wild values here
    p[p.<=-3] .= [-3]
    p[p.>4] .= [3]
    
    for j in 1:nitemsgrm
        a[j, 1] = p[(j - 1) * nParmsPerItemGRM + 1]
        a_z[j, 1] = p[nitems * nParmsPerItemGRM + (j - 1) * nParmsPerItem2PL + 1]
        b_z[j, 1] = p[nitems * nParmsPerItemGRM + (j - 1) * nParmsPerItem2PL + 2]
        for k in 1:(ncatgrm-1)
            b[j, k] = p[(j - 1) * nParmsPerItemGRM + 1 + k]
        end
    end

    rho_exp = p[nitems * nParmsPerItemGRM + nitems * nParmsPerItem2PL + 1]
    rho = exp(rho_exp) / (1 + exp(rho_exp))
    itemtrace = trace_line_pts(a, b, a_z, b_z, theta);
    expected = zeros(size(r, 1))

    ## Create all of the posterior estimates here
    posterior = zeros(size(theta)[1])
    for i in 1:size(theta)[1]
        posterior[i] = multivariate_normal_pdf(theta[i,:], [0,0], [1.0 rho; rho 1.0])
    end
    posterior_orig = posterior  
    for i in 1:size(r)[1]
        posterior = posterior_orig
        for item in 1:size(testdata)[2]
            x = Int(r[i, item])   
            posterior = posterior .* itemtrace[item, :,x+1] 
        end
        expected[i] = sum(posterior)
    end
    #@info "params are: $p"
    ## Check for 0 values in the expected values
    expected[expected.==0] .= [1e-10]
    expected[expected.<0] .= [1e-10]    
    l = -1 * sum( r[:,15] .* log.(expected))
    if isnan(l)
        @info "params are: $p"
        error("LL is NaN")
    end
    if isinf(l)
        @info "params are $p"
        error("LL is inf")
    end
    #@info "LL val is: $l"
    return(l)
    
end


## Now call the ll model
ll_grm_ip(p, testitems, theta)

@time ll_grm_ip(p, testitems, theta)


## Now try to optimize this
@time h = optimize(z -> ll_grm_ip(z, testitems, theta),p,LBFGS(),Optim.Options(g_tol = 1e-3, iterations=350_000, show_trace=true, show_every=5)) ## This took 21479 seconds

h_min = Optim.minimizer(h);
println(h_min)
h_min = DataFrame(h_min)

@time j = optimize(z -> ll_grm_ip(z, testitems, theta),p,Optim.Options(g_tol = 1e-3, iterations=350_000, show_trace=true, show_every=100))
p2 = Optim.minimizer(j);
@time j2 = optimize(z -> ll_grm_ip(z, testitems, theta),p2,Optim.Options(g_tol = 1e-3, iterations=25_000, show_trace=true, show_every=100))


## Now try NM with range restiriction
lowerBound = p .- [-6]
upperBound = p .+  [6]

## Now obtain the fscores
Dep = Optim.minimizer(j);
Dep = reshape(Dep, length(p), 1)

a_grm = Dep[1:nitems*nParmsPerItemGRM:nParmsPerItemGRM] # Select a parameters from GRM
a_grm = DataFrame(a_grm, nrow=nitems)

b_grm = zeros(nitems, ncatgrm - 1) # Select b parameters from GRM
for i in 1:nitems
    b_grm[i, :] .= Dep[(ncatgrm*i - 1):(ncatgrm*i)]
end
b_grm = DataFrame(b_grm, nrow=nitems, ncol=ncatgrm - 1)

a_2PL = Dep[nitems*nParmsPerItemGRM + 1:nitems*nParmsPerItemGRM + nitems*nParmsPerItem2PL:2] # Select a parameters from 2PL
a_2PL = DataFrame(a_2PL, nrow=nitems)

b_2PL = Dep[nitems*nParmsPerItemGRM + 2:nitems*nParmsPerItemGRM + nitems*nParmsPerItem2PL:2] # Select b parameters from 2PL
b_2PL = DataFrame(b_2PL, nrow=nitems)

rho = exp(Dep[end]) / (1 + exp(Dep[end]))

itemtrace2PL = trace_line_pts_2PL(a_2PL, b_2PL, theta)
itemtraceGRM = trace_line_pts_grm(a_grm, b_grm, theta)
itemtrace = trace_line_pts(a_grm, b_grm, a_2PL, b_2PL, theta)

qpoints = theta
prior = mvnormal(qpoints, [0, 0], [1 rho; rho 1])

function score(pattern)
    lhood = ones(length(qpoints.theta1))
    for item in 1:nitems
        if pattern[item] == 0
            lhood .*= itemtrace[1][item, :]
        elseif pattern[item] == 1
            lhood .*= itemtrace[2][item, :]
        elseif pattern[item] == 2
            lhood .*= itemtrace[3][item, :]
        elseif pattern[item] == 3
            lhood .*= itemtrace[4][item, :]
        end
    end
    return lhood
end

for respondent in 1:size(data_out, 1)
    pattern = response[respondent, findall(x -> occursin(r"x\d+", x), names(response))]
    lhood = score(pattern)
    
    eap2PL_Hurdle = sum(lhood .* prior .* qpoints.theta1) / sum(lhood .* prior)
    se2PL_Hurdle = sqrt(sum(lhood .* prior .* (qpoints.theta1 .- eap2PL_Hurdle).^2) / sum(lhood .* prior))
    data_out.eap2PL_Hurdle[respondent] = eap2PL_Hurdle
    data_out.se2PL_Hurdle[respondent] = se2PL_Hurdle
    
    eapGRM_Hurdle = sum(lhood .* prior .* qpoints.theta2) / sum(lhood .* prior)
    seGRM_Hurdle = sqrt(sum(lhood .* prior .* (qpoints.theta2 .- eapGRM_Hurdle).^2) / sum(lhood .* prior))
    data_out.eapGRM_Hurdle[respondent] = eapGRM_Hurdle
    data_out.seGRM_Hurdle[respondent] = seGRM_Hurdle
end

