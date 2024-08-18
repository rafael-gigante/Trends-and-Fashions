module plots

using Plots
using LaTeXStrings
using DataFrames
using GLM
theme(:default)

"""
    plot_measures(iteration, entropy, percolation, suscept, N, L, β,  n, l, b)

Plot the entropy, percolation and susceptibility given a value of 'n ∈ N', 'l ∈ L' and 'b ∈ β' that has been used in the simulation.

"""
function plot_measures(iteration, entropy, percolation, suscept, N, L, β, ensembles, n, l, b)
    n_index = findfirst(==(n), N)
    l_index = findfirst(==(l), L)
    b_index = findfirst(==(b), β)

    pl = plot(iteration[n_index][l_index][b_index][1:ensembles], entropy[n_index][l_index][b_index][1:ensembles], xlabel="time [iterations]", ylabel=L"$S$", label=false)
    hline!([log10(l)], line=:dash, linewidth=1, color=:black, label=false)
    plot!(minorgrid=true, title="N=$n, L=$l, β =$b", framestyle=:box, legend=:right )
    ylims!(0, log10(l) + 0.02)
    display(pl)

    pl = plot(iteration[n_index][l_index][b_index][1:ensembles], percolation[n_index][l_index][b_index][1:ensembles], xlabel="time [iterations]", ylabel=L"$P_c/N$", label=false)
    hline!([1], line=:dash, linewidth=1, color=:black, label=false)
    plot!(minorgrid=true, title="N=$n, L=$l, β =$b", framestyle=:box, legend=:right )
    ylims!(0, 1.02)
    display(pl)

    pl = plot(iteration[n_index][l_index][b_index][1:ensembles], suscept[n_index][l_index][b_index][1:ensembles], xlabel="time [iterations]", ylabel=L"$S_c/N^2$", label=false)
    plot!(minorgrid=true, title="N=$n, L=$l, β =$b", framestyle=:box, legend=:right )
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
    plot_τ_fixed_l(τ, N, L, β, l)

Make a scatter plot of log(τ) vs log(N) for each β with a fixed value of l ∈ L and fit the data with a linear regression.

"""
function plot_τ_fixed_l(τ, N, L, β, l)
    l_index = findfirst(==(l), L)

    pl = plot(xlabel=L"$log(N)$", ylabel=L"$log(\tau)$", title="L=$l", minorgrid=true, framestyle=:box)
    for b in eachindex(β)
        temp_τ = []
        for n in eachindex(N)
            push!(temp_τ, τ[n][l_index][b])
        end
        data, slope, intercept = fit_data(log10.(N), log10.(temp_τ))

        scatter!(data.X, data.Y, label=L"$\beta = %$(β[b])$", color=palette(:default)[b])
        plot!(data.X, data.Y_pred, label=L"$Y = %$intercept + %$slope * X$", line=:dash, color=palette(:default)[b])
    end
    display(pl)
end

"""
    plot_τ_fixed_β(τ, N, L, β, b)

Make a scatter plot of log(τ) vs log(N) for each L with a fixed value of b ∈ β and fit the data with a linear regression.

"""
function plot_τ_fixed_β(τ, N, L, β, b)
    b_index = findfirst(==(b), β)

    pl = plot(xlabel=L"$log(N)$", ylabel=L"$log(\tau)$", title="β=$b", minorgrid=true, framestyle=:box)
    for l in eachindex(L)
        temp_τ = []
        for n in eachindex(N)
            push!(temp_τ, τ[n][l][b_index])
        end
        data, slope, intercept = fit_data(log10.(N), log10.(temp_τ))

        scatter!(data.X, data.Y, label=L"$L = %$(L[l])$", color=palette(:default)[l])
        plot!(data.X, data.Y_pred, label=L"$Y = %$intercept + %$slope * X$", line=:dash, color=palette(:default)[l])
    end
    display(pl)
end

"""
    fermi_function(x, β)

Returns the value of the Fermi-Dirac function given a value 'x' and a temperature 'β'.

"""
function fermi_function(x, β)
    return  1 / (1 + exp((-β*x)))
end

"""
    plot_fermi_function(β)

Plot the Fermi-Dirac function for x ∈[-1,1] for a given set of 'β'.

"""
function plot_fermi_function(β)
    x = -1:0.001:1
    pl = plot()
    for j in eachindex(β)
        y = []
        for i in x
            push!(y, fermi_function(i, β[j]))
        end
        plot!(x, y, label="β = $(β[j])")
    end
    plot!(xlabel=L"x", ylabel=L"Probabilidade $q_{ij}$", minorgrid=true, framestyle=:box,)
    ylims!(-0.02, 1.02)
    xlims!(-1,1)
    display(pl)
end

"""
    plot_labels_evolution(iteration, labels_evolution, N, L, β, ensembles, n_greatest, n, l, b)

Plot the evolution of the 'n_greatest' labels of a random ensamble for given a value of 'n ∈ N', 'l ∈ L' and 'b ∈ β' that has been used in the simulation.

"""
function plot_labels_evolution(iteration, labels_evolution, N, L, β, ensembles, n_greatest, n, l, b)
    n_index = findfirst(==(n), N)
    l_index = findfirst(==(l), L)
    b_index = findfirst(==(b), β)

    pl = plot(xlabel="time [iterations]", ylabel=L"$P_c/N$", title="N=$n, L=$l, β =$b", minorgrid=true, framestyle=:box)
    a = rand(1:ensembles)
    for i in 1:n_greatest
        plot!(iteration[n_index][l_index][b_index][a], labels_evolution[n_index][l_index][b_index][a][i], label="$(i)° tendência")
    end
    hline!([1], line=:dash, linewidth=1, color=:black, label=false)
    ylims!(0, 1.02)
    display(pl)
end

"""
    plot_τ_fixed_β(τ, N, L, β, l)

Make a scatter plot of log(τ) vs log(N) for each L with a fixed value of b ∈ β and fit the data with a linear regression.

"""
function plot_τ_fixed_β(τ, N, L, β, b)
    b_index = findfirst(==(b), β)

    pl = plot(xlabel=L"$log(N)$", ylabel=L"$log(\tau)$", title="β=$b", minorgrid=true, framestyle=:box)
    for l in eachindex(L)
        temp_τ = []
        for n in eachindex(N)
            push!(temp_τ, τ[n][l][b_index])
        end
        data, slope, intercept = fit_data(log10.(N), log10.(temp_τ))

        scatter!(data.X, data.Y, label=L"$L = %$(L[l])$", color=palette(:default)[l])
        plot!(data.X, data.Y_pred, label=L"$Y = %$intercept + %$slope * X$", line=:dash, color=palette(:default)[l])
    end
    display(pl)
end

"""
    plot_τ_over_N(τ, N, L, β, l)

Make a plot of τ/N vs N for a given value of l ∈ L.

"""
function plot_τ_over_N(τ, N, L, β, l)
    l_index = findfirst(==(l), L)

    pl = plot(xlabel=L"N", ylabel=L"$\tau/N$", minorgrid=true, framestyle=:box)
    for b in eachindex(β)
        temp_τ = []
        for n in eachindex(N)
            push!(temp_τ, τ[n][l_index][b]/N[n])
        end

        plot!(N, temp_τ, label=L"$β = %$(β[b])$", color=palette(:default)[b])
    end
    display(pl)
end

end