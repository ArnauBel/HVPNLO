
function jackknife_resampling(data)
    N = length(data)
    means = zeros(N)
    
    for i in 1:N
        sample = deleteat!(copy(data), i)  # Leave-one-out sample
        means[i] = mean(sample)
    end
    
    mean_jack = mean(means)
    std_jack = sqrt((N - 1) * mean((means .- mean_jack) .^ 2))
    
    return [mean_jack, std_jack]
end

function bootstrap_resampling(data, B=1000)
    N = length(data)
    means = zeros(B)
    
    for i in 1:B
        resample = data[rand(1:N, N)]  # Resample with replacement
        means[i] = mean(resample)
    end
    
    mean_boot = mean(means)
    std_boot = std(means)  # Standard deviation of bootstrap samples
    
    return [mean_boot, std_boot]
end

function jackknife_err(jakknife_sample,∆dof=0)
    N    = length(jakknife_sample)
    Mean = mean(jakknife_sample)
    Err  = sqrt(N*(N-1)/(N-∆dof)) * sqrt(mean((jakknife_sample.-Mean).^2))
    return [Mean, Err]
end