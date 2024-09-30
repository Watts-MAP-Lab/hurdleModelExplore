using Distributions, StatsBase, DataFrames


## I am going to try to take this function and write it myself so it actually works
function simIrtMod2PL(a::Float16, b::Float16, N::Int, θ_range::Tuple{Float16, Float16}=(-3.0, 3.0))
    ## First simulate θ vals
    Θ = rand(Normal(0.0, 1.0), N) .* (θ_range[2] - θ_range[1]) .+ θ_range[1]; ## Note \Theta; not \theta!
    
    ## Now create output responses
    responses = zeros(N, size(a)[1]);
    ## Now sample through item and participant
    for i in 1:N
        for j in 1:size(a)[1]
            p = 1 / (1 + exp(-a[j] * (Θ[i] - b[j])))
            responses[i, j] = rand(Bernoulli(p))
        end
    end
    ## Now return the output
    df = DataFrame(Θ = Θ, responses = responses)
    return(df)
end

## Now build a multivivariate 2PL model?
function simIrtModMV2PL(μ, Σ, a, b, N)
    ## first create the Θ values
    Θ = rand(MvNormal(μ, Σ), N)
    ## Now get the response probability
    ## DO this in a loop
    responses = zeros(N, size(a)[2])
    for i in 1:N
        for j in 1:size(a)[2]
            p = sum(a[:,j] .* Θ[:,i]) - b[j]
            p = 1 / (1 + exp(-p))
            responses[i,j] = rand(Bernoulli(p))
        end
    end

    ## Now prep the output
    df = DataFrame(Θ = Θ, responses = responses)
    return(df)
end


## 2PL prob here
function item_response_probability(a, b, θ)
    #a = a[j]
    #b = b[j]
    return 1 / (1 + exp(-a * (θ - b)))
end

## Now do graded response model here
function response_probability_grm(θ, a, b, k)
    # θ: latent trait of the respondent
    # a: discrimination parameter of the item
    # b: vector of thresholds (cut points) for the item
    # k: the category for which the probability is being computed (k = 1, 2, ..., K)
    
    if k == 1
        # For the lowest category, P(Y = 1) = 1 - P(Y ≥ 2)
        P_k = 1 - logistic(a * (θ - b[1]))
    elseif k == length(b) + 1
        # For the highest category, P(Y = K) = P(Y ≥ K)
        P_k = logistic(a * (θ - b[end]))
    else
        # For intermediate categories, P(Y = k) = P(Y ≥ k) - P(Y ≥ k + 1)
        P_k = logistic(a * (θ - b[k - 1])) - logistic(a * (θ - b[k]))
    end
    
    return P_k
end


function generate_item_parameters2PL(num_items, num_categories)
    a = rand(Uniform(0.5, 2.0), num_items)  # discrimination parameters
    b = rand(Uniform(-2.0, 2.0), num_items)  # thresholds for each item
    return a, b
end


function generate_item_parametersGRM(num_items, num_categories)
    a = rand(Uniform(0.5, 2.0), num_items)  # discrimination parameters
    b = [sort(rand(Uniform(-2.0, 2.0), num_categories - 1)) for _ in 1:num_items]  # thresholds for each item
    return a, b
end

function logistic(x)
    return 1.0 / (1.0 + exp(-x))
end

function simulate_hurdle_responses(a, b, a_z, b_z, μ, Σ, N,k)
    ## First create the Θ values
    Θ = rand(MvNormal(μ, Σ), N);

    ## Data need to be generated from a mixture model which takes the following form:
    # T_{j} = (Y_{ij}=k|Θ_{0i},Θ_{1i},α_j)=[T_{0j}=(Y_{ij}=0|Θ_{0i},α_j)] X
    #       [T_{0j}(Y_{ij}≠0|Θ_{0i},α_j) T_1j]

    ## Now grab the 2PL probailities first -- this is the symptom susceptibility portion of the model
    responses_sus = zeros(N, length(a));
    for i in 1:N
        for j in 1:length(a)
            responses_sus[i, j] = item_response_probability(a[j], b[j], Θ[1,i]) ## This gives the probability of an endorsment
        end
    end
    ## Now create the inverse of this -- the probability of no endorsment
    responses_neg = 1 .- responses_sus;
    ## Now I need to obtain the probability of responding to each positive category
    responses_cats = zeros(N, length(a), k-1);
    ## Now loop through every participant, item, and response category
    for i in 1:N
        for j in 1:length(a)
            ## Give the probability of an endorsment
            prob_endorse = 1 - response_probability(a[j], b[j], θ[1,i])    
            for l in 1:k-1
                response_cats[i,j,l] = response_probability_grm(Θ[2,i], a_z[j], b_z[j,:], l) * prob_endorse 
            end
        end
    end

end