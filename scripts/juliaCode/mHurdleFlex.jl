using Combinatorics, Distributions, CSV, DataFrames, Optim, LinearAlgebra, Suppressor

## First load the data
in_resp = ARGS[1]
in_tabs = ARGS[2]
#in_resp = "/home/arosen/Documents/hurdleModelExplore/data/testOSF.csv"
#in_tabs = "/home/arosen/Documents/hurdleModelExplore/data/testOSF2.csv"
data_out = CSV.File(in_resp) |> DataFrame
data_out2 = CSV.File(in_tabs) |> DataFrame

## Now obtain all of the names in the data frame
println("Processing: ", in_resp)


## Create prob estimate for GRM portion
function trace_line_pts_grm(a, b, theta)
    n_categories = size(b)[2] + 1
    itemtraceGRM = zeros(size(a)[1], size(theta)[1], n_categories)
    for j in 1:size(a, 1)
        for k in 1:n_categories
            if k == 1
                itemtraceGRM[j, :,1] .= 1 .- exp.(a[j, :] .* (theta[:, 2] .- b[j, k])) ./ (1 .+ exp.(a[j, :] .* (theta[:, 2] .- b[j, k])))
            elseif k > 1 && k < n_categories 
                itemtraceGRM[j, :,k] .= exp.(a[j, :] .* (theta[:, 2] .- b[j, k-1])) ./ (1 .+ exp.(a[j, :] .* (theta[:, 2] .- b[j, k-1]))) .- exp.(a[j, :] .* (theta[:, 2] .- b[j, k])) ./ (1 .+ exp.(a[j, :] .* (theta[:, 2] .- b[j, k])))
            elseif k == n_categories
                itemtraceGRM[j, :,k] .= exp.(a[j, :] .* (theta[:, 2] .- b[j, k-1])) ./ (1 .+ exp.(a[j, :] .* (theta[:, 2] .- b[j, k-1])))
            end
        end
    end
    return itemtraceGRM
end

## Create prob estimate for 2PL portion
function trace_line_pts_2PL(a_z, b_z, theta)
    nitems = size(a_z)[1];
    itemtrace2PL = zeros(nitems, size(theta)[1]);
    for j in 1:nitems
        itemtrace2PL[j, :] .= 1 .- exp.(a_z[j, :] .* (theta[:, 1] .- b_z[j, 1])) ./ (1 .+ exp.(a_z[j, :] .* (theta[:, 1] .- b_z[j, 1])))
    end
    return itemtrace2PL
end

## Combine 2PL and GRM model estimates here
function trace_line_pts(a, b, a_z, b_z, theta)
    n_categories = size(b)[2]+2
    itemtrace = zeros(size(a)[1],size(theta)[1],n_categories)
    for j in 1:size(a)[1]
        for k in 0:n_categories-1
            #print(k)
            if k == 0
                itemtrace[j, :,k+1] .= 1 .- trace_line_pts_2PL(a_z, b_z, theta)[j, :]
            elseif k > 0
                itemtrace[j,:,k+1] .= trace_line_pts_2PL(a_z, b_z, theta)[j, :] .* trace_line_pts_grm(a, b, theta)[j, :,k]
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


function ll_grm_ip(p, testdata, theta, r)
    
    nParamsPerItemGRM = maximum(maximum(eachcol(testdata))); 
    nParamsPerItem2PL = 2
    ncatgrm = maximum(maximum(eachcol(testdata))) 
    nitemsgrm=nitems=size(testdata)[2]
    
    a = fill(-1.0, nitems, 1)
    b = fill(-1.0, nitems, nParamsPerItemGRM-1)
    a_z = fill(-1.0, nitems, 1)
    b_z = fill(-1.0, nitems, 1)
    rho_exp = 0

    ## Fix wild values here
    p[p.<=-8] .= [-8]
    p[p.>8] .= [8]
    
    ## Now get the flag index here
    flag_index = fill(0, size(p)[1], 1);

    for j in 1:nitemsgrm
        a[j, 1] = p[(j - 1) * nParmsPerItemGRM + 1]
        flag_index[(j - 1) * nParmsPerItemGRM + 1] = flag_index[(j - 1) * nParmsPerItemGRM + 1] + 1 
        a_z[j, 1] = p[nitems * nParmsPerItemGRM + (j - 1) * nParmsPerItem2PL[1]+ 1]
        flag_index[nitems * nParmsPerItemGRM + (j - 1) * nParmsPerItem2PL[1] + 1] = flag_index[nitems * nParmsPerItemGRM + (j - 1) * nParmsPerItem2PL[1] + 1] + 1
        b_z[j, 1] = p[nitems * nParmsPerItemGRM + (j - 1) * nParmsPerItem2PL[1] + 2]
        flag_index[nitems * nParmsPerItemGRM + (j - 1) * nParmsPerItem2PL[1] + 2] = flag_index[nitems * nParmsPerItemGRM + (j - 1) * nParmsPerItem2PL[1] + 2] + 1
        for k in 1:(ncatgrm-1)
            b[j, k] = p[(j - 1) * nParmsPerItemGRM + 1 + k]
            flag_index[(j - 1) * nParmsPerItemGRM + 1 + k] = flag_index[(j - 1) * nParmsPerItemGRM + 1 + k] + 1
        end
    end

    rho_exp = p[nitems * nParmsPerItemGRM + nitems * nParamsPerItem2PL + 1]
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
        for item in 1:nitems
            x = Int(r[i, item])   
            posterior = posterior .* itemtrace[item, :,x+1] 
        end
        expected[i] = sum(posterior)
    end
    #@info "params are: $p"
    ## Check for 0 values in the expected values
    expected[expected.==0] .= [1e-10]
    expected[expected.<0] .= [1e-10]    
    l = -1 * sum( r[:,end] .* log.(expected))
    if isnan(l)
        @info "params are: $p"
        error("LL is NaN")
    end
    if isinf(l)
        @info "params are $p"
        error("LL is inf")
    end
    #@info "LL val is: $l"
    #@info "params are $p"

    return(l)
    
end

## Basic shop keeping here
nitems = size(data_out)[2];
nParmsPerItemGRM = maximum(maximum(eachcol(data_out))); ## Finds the maximum number of response options and then subtract by one for zero portion
nParmsPerItem2PL = [2]; ## By definition, only 2 params in 2 PL model

## Now create all theta points for Quadrature integration nonsense
theta1 = -6:0.25:6  # Quadrature points for theta1
theta2 = -6:0.25:6
theta = collect(Combinatorics.permutations(theta1,2));
d = collect(Iterators.product(theta1, theta2));
d = reshape(d, 2401);
theta = [d[i][j] for i in 1:length(d), j in 1:length(d[1])];

## Now initialize random start points for GRM postion of model
a = rand(nitems);
b = rand(nitems, nParmsPerItemGRM-1);

## Now intitialize random start points for the 2PL portion of model
a_z = rand(size(data_out)[2]);
b_z = rand(size(data_out)[2]);

rho_exp = 0.2;

paramGRM = nitems*nParmsPerItemGRM
param2PL = nitems*nParmsPerItem2PL
totalP = paramGRM .+ param2PL
totalP = totalP[1] + 1
p = zeros(totalP) 
nitemsgrm = nitems
ncatgrm = maximum(maximum(eachcol(data_out)))
for j in 1:nitems 
  p[(j-1)*nParmsPerItemGRM + 1] = a[j]
  p[nitemsgrm[1]*nParmsPerItemGRM[1] + (j-1)*nParmsPerItem2PL[1] + 1] = a_z[j]
  p[nitemsgrm[1]*nParmsPerItemGRM[1] + (j-1)*nParmsPerItem2PL[1] + 2] = b_z[j]
  for k in 1:(ncatgrm-1)
    p[(j-1)*nParmsPerItemGRM[1] + 1 + k] = b[j,k]
  end 
end
p[nitems[1]*nParmsPerItemGRM[1] + nitems[1]*nParmsPerItem2PL[1] + 1] = rho_exp

## Now optimize these values

## First run single shot ll test
println(ll_grm_ip(p, data_out, theta, data_out2))

## Now optimize
## Now write the csv here
#outputRand = rand(Int)
#outputFile = ("/tmp/")

@time h = optimize(z -> ll_grm_ip(z, data_out, theta, data_out2),p,LBFGS(),Optim.Options(g_tol = 1e-3, iterations=350_000, show_trace=true, show_every=5))

println(Optim.minimizer(h))
