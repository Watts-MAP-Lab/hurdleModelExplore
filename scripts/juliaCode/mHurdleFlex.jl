using Combinatorics, Distributions, CSV, DataFrames, Optim, LinearAlgebra#, Suppressor

## First load the data
#in_resp = ARGS[1]
in_tabs = ARGS
#in_resp = "/home/arosen/Documents/hurdleModelExplore/data/testOSF.csv"
#in_tabs = "/home/arosen/Documents/hurdleModelExplore/data/testOSF2.csv"
#data_out = CSV.File(in_resp) |> DataFrame
data_out2 = CSV.File(in_tabs) |> DataFrame

## Now obtain all of the names in the data frame
println("Processing: ", in_tabs)


## Create prob estimate for GRM portion
function trace_line_pts_grm(a, b, theta)
    n_items, n_theta = size(a, 1), size(theta, 1)
    n_categories = size(b, 2) + 1
    itemtraceGRM = zeros(n_items, n_theta, n_categories)
    theta_col2 = theta[:, 2]  # Extract the second column of theta once

    for j in 1:n_items
        for k in 1:n_categories
            if k == 1
                diff_b1 = a[j, :] .* (theta_col2 .- b[j, k])
                exp_b1 = exp.(diff_b1)
                itemtraceGRM[j, :, 1] .= 1 .- exp_b1 ./ (1 .+ exp_b1)
            elseif k < n_categories
                diff_bk1 = a[j, :] .* (theta_col2 .- b[j, k-1])
                diff_bk = a[j, :] .* (theta_col2 .- b[j, k])
                exp_bk1 = exp.(diff_bk1)
                exp_bk = exp.(diff_bk)
                itemtraceGRM[j, :, k] .= exp_bk1 ./ (1 .+ exp_bk1) .- exp_bk ./ (1 .+ exp_bk)
            else
                diff_bk1 = a[j, :] .* (theta_col2 .- b[j, k-1])
                exp_bk1 = exp.(diff_bk1)
                itemtraceGRM[j, :, k] .= exp_bk1 ./ (1 .+ exp_bk1)
            end
        end
    end

    return itemtraceGRM
end
## Create prob estimate for 2PL portion
function trace_line_pts_2PL(a_z, b_z, theta)
    nitems = size(a_z, 1)
    ntheta = size(theta, 1)
    itemtrace2PL = zeros(nitems, ntheta, 2)
    theta_col1 = theta[:, 1]  # Extract the first column of theta once

    for j in 1:nitems
        diff_bz = a_z[j, :] .* (theta_col1 .- b_z[j, 1])
        exp_bz = exp.(diff_bz)
        denom = 1 .+ exp_bz

        itemtrace2PL[j, :, 1] .= 1 .- exp_bz ./ denom  # Probability of "no"
        itemtrace2PL[j, :, 2] .= exp_bz ./ denom      # Probability of "yes"
    end

    return itemtrace2PL
end


## Combine 2PL and GRM model estimates here
function trace_line_pts(a, b, a_z, b_z, theta)
    n_categories = size(b, 2) + 2
    n_items = size(a, 1)
    n_theta = size(theta, 1)

    itemtrace = zeros(n_items, n_theta, n_categories)

    # Precompute the results from trace_line_pts_2PL and trace_line_pts_grm
    itemtrace2PL = trace_line_pts_2PL(a_z, b_z, theta)
    itemtraceGRM = trace_line_pts_grm(a, b, theta)

    for j in 1:n_items
        for k in 0:n_categories-1
            if k == 0
                # Use the precomputed 2PL results for "no" probability
                itemtrace[j, :, k+1] .= itemtrace2PL[j, :, 1]
            else
                # Use the precomputed 2PL "yes" and GRM results
                itemtrace[j, :, k+1] .= itemtrace2PL[j, :, 2] .* itemtraceGRM[j, :, k]
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
        println(Σ)
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

function ll_grm_ip(p, data_in, theta)
    # First remove the count column
    testdata = data_in[:, 1:end-1]
    r = data_in[:, end]
    
    nParamsPerItemGRM = maximum(maximum(eachcol(testdata)))
    nParamsPerItem2PL = 2
    ncatgrm = maximum(maximum(eachcol(testdata)))
    nitemsgrm = nitems = size(testdata, 2)
    
    # Initialize parameters a, b, a_z, b_z
    a = fill(-1.0, nitems, 1)
    b = fill(-1.0, nitems, nParamsPerItemGRM - 1)
    a_z = fill(-1.0, nitems, 1)
    b_z = fill(-1.0, nitems, 1)
    
    # Extract parameters from p
    for j in 1:nitemsgrm
        idx_grm = (j - 1) * nParamsPerItemGRM + 1
        idx_2pl = nitems * nParamsPerItemGRM + (j - 1) * nParamsPerItem2PL
        a[j, 1] = p[idx_grm]
        a_z[j, 1] = p[idx_2pl + 1]
        b_z[j, 1] = p[idx_2pl + 2]
        
        for k in 1:(ncatgrm - 1)
            b[j, k] = p[idx_grm + k]
        end
    end
    
    # Sort difficulty parameters for each item
    b = sort(b, dims=2)
    
    # Calculate rho
    rho_exp = p[nitems * nParamsPerItemGRM + nitems * nParamsPerItem2PL + 1]
    rho = exp(rho_exp) / (1 + exp(rho_exp))
    rho = clamp(rho, 1e-10, 0.99999999)  # Avoid edge cases for rho
    
    # Ensure all a & a_z are positive
    a = abs.(a)
    a_z = abs.(a_z)
    
    # Precompute itemtrace
    itemtrace = trace_line_pts(a, b, a_z, b_z, theta)
    
    # Precompute posterior distribution using multivariate normal pdf
    posterior_orig = [multivariate_normal_pdf(theta[i, :], [0, 0], [1.0 rho; rho 1.0]) for i in 1:size(theta, 1)]
    
    # Initialize expected value array
    expected = zeros(size(r))
        for i in 1:size(r)[1]
        posterior = posterior_orig
        for item in 1:nitems
            x = Int(testdata[i, item])   
            posterior = posterior .* itemtrace[item, :,x+1] 
        end
        expected[i] = sum(posterior)
    end
    
    # Compute log likelihood, adding a small constant to avoid log(0)
    expected = max.(expected, 1e-10)
    l = -sum(r .* log.(expected))

    return l
end


## Basic shop keeping here
nitems = size(data_out2)[2]-1;
nParmsPerItemGRM = maximum(maximum(eachcol(data_out2[:,1:end-1]))); ## Finds the maximum number of response options and then subtract by one for zero portion
nParmsPerItem2PL = [2]; ## By definition, only 2 params in 2 PL model
## Now create all theta points for Quadrature integration nonsense
theta1 = -6:0.25:6  # Quadrature points for theta1
theta2 = -6:0.25:6
theta = collect(Combinatorics.permutations(theta1,2));
d = collect(Iterators.product(theta1, theta2));
d = reshape(d, 2401);
theta = [d[i][j] for i in 1:length(d), j in 1:length(d[1])];

## Now initialize random start points for GRM postion of model
a = fill(2, nitems);
b = rand(nitems, nParmsPerItemGRM-1);

## Now intitialize random start points for the 2PL portion of model
a_z = fill(2, nitems);
b_z = fill(-1, nitems);

rho_exp = 0.2;

paramGRM = nitems*nParmsPerItemGRM
param2PL = nitems*nParmsPerItem2PL
totalP = paramGRM .+ param2PL
totalP = totalP[1] + 1
p = zeros(totalP) 
nitemsgrm = nitems
ncatgrm = maximum(maximum(eachcol(data_out2[:,1:end-1])))
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
println(ll_grm_ip(p, data_out2, theta))
@time ll_grm_ip(p, data_out2, theta)

## Now optimize
h = optimize(z -> ll_grm_ip(z, data_out2, theta),p,LBFGS(),Optim.Options(g_tol = 1e-3, iterations=350_000, show_trace=true, show_every=5))

println(Optim.minimizer(h))