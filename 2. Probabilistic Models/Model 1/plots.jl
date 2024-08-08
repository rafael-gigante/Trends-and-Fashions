module plots

using Plots
using LaTeXStrings
using DataFrames
using GLM
theme(:default)

"""
    plot_measures(iteration, entropy, percolation, suscept, N, L, ensembles, n, l)

Plot the entropy, percolation and susceptibility given a value of 'n ∈ N' and 'l ∈ L' that has been used in the simulation.

"""
function plot_measures(iteration, entropy, percolation, suscept, N, L, ensembles, n, l)
    n_index = findfirst(==(n), N)
    l_index = findfirst(==(l), L)

    pl = plot(iteration[n_index][l_index][1:ensembles], entropy[n_index][l_index][1:ensembles], xlabel="time [iterations]", ylabel=L"$S$", label=false)
    hline!([log10(l)], line=:dash, linewidth=1, color=:black, label=false)
    plot!(minorgrid=true, title="N=$n, L=$l", framestyle=:box, legend=:right )
    ylims!(0, log10(l) + 0.02)
    display(pl)

    pl = plot(iteration[n_index][l_index][1:ensembles], percolation[n_index][l_index][1:ensembles], xlabel="time [iterations]", ylabel=L"$P_c/N$", label=false)
    hline!([1], line=:dash, linewidth=1, color=:black, label=false)
    plot!(minorgrid=true, title="N=$n, L=$l", framestyle=:box, legend=:right )
    ylims!(0, 1.02)
    display(pl)

    pl = plot(iteration[n_index][l_index][1:ensembles], suscept[n_index][l_index][1:ensembles], xlabel="time [iterations]", ylabel=L"$S_c/N^2$", label=false)
    plot!(minorgrid=true, title="N=$n, L=$l", framestyle=:box, legend=:right )
    display(pl)
end

"""
    fit_data(x, y)

Receive two arrays and fits the data with a linear regression.

Returns:
- 'slope': Float with the value of the angular coefficient
- 'intercept': Float with the value of the linear coefficient
"""

function fit_data(x, y)
    # Create DataFrame
    data = DataFrame(
        X = x,
        Y  = y
    )
    # Fit the model
    model = lm(@formula(Y ~ X), data)
    # Get the fitted values
    data.Y_pred = predict(model)
    coefficients = coef(model)
    intercept = round(coefficients[1], digits=2)
    slope = round(coefficients[2], digits=2)
    return data, slope, intercept
end

"""
    plot_τ_fixed_l(τ, N, L)

Make a scatter plot of log(τ) vs log(N) for each value of l ∈ L and fit the data with a linear regression.

"""
function plot_τ_fixed_L(τ, N, L)

    pl = plot(xlabel=L"$log(N)$", ylabel=L"$log(\tau)$", minorgrid=true, framestyle=:box)
    for l in eachindex(L)
        temp_τ = []
        for n in eachindex(N)
            push!(temp_τ, τ[n][l])
        end
        data, slope, intercept = fit_data(log10.(N), log10.(temp_τ))

        scatter!(data.X, data.Y, label=L"$L = %$(L[l])$", color=palette(:default)[l])
        plot!(data.X, data.Y_pred, label=L"$Y = %$intercept + %$slope * X$", line=:dash, color=palette(:default)[l])
    end
    display(pl)
end

"""
    plot_τ_fixed_N(τ, N, L)

Make a scatter plot of log(τ) vs log(L) for each value of n ∈ N and fit the data with a linear regression.

"""
function plot_τ_fixed_N(τ, N, L)

    pl = plot(xlabel=L"$log(L)$", ylabel=L"$log(\tau)$", minorgrid=true, framestyle=:box)
    for n in eachindex(N)
        temp_τ = []
        for l in eachindex(L)
            push!(temp_τ, τ[n][l])
        end
        data, slope, intercept = fit_data(log10.(L), log10.(temp_τ))

        scatter!(data.X, data.Y, label=L"$N = %$(N[n])$", color=palette(:default)[n])
        plot!(data.X, data.Y_pred, label=L"$Y = %$intercept + %$slope * X$", line=:dash, color=palette(:default)[n])
    end
    display(pl)
end

"""
    plot_τ_NL(τ, N, L)

Make a scatter plot of log(τ) vs log(N)log(L) and fit the data with a linear regression.

"""
function plot_τ_NL(τ, N, L)

    pl = plot(minorgrid=true, framestyle=:box, xlabel=L"log(N)log(L)", ylabel=L"log($\tau$)")

    x_axis = []
    y_axis = []
    for l in eachindex(L)
        x_temp = []
        y_temp = []
        for n in eachindex(N)
            push!(x_temp, log10.(N[n])*log10.(L[l]))
            push!(y_temp, log10.(τ[n][l]))
        end
        append!(x_axis, x_temp)
        append!(y_axis, y_temp)
        scatter!(x_temp, y_temp, label=L"log(L) = %$(log10.(L[l]))")
    end 

    display(pl)

end

end