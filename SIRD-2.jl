### A Pluto.jl notebook ###
# v0.14.8

using Markdown
using InteractiveUtils

# ╔═╡ 9f18ccd0-8e25-11eb-20e6-ab941f663550
begin
	import Pkg
    Pkg.activate(".")
    Pkg.add("PlutoUI")
	Pkg.add("Plots")
	Pkg.add("Random")
	Pkg.add("Distributions")
	
	using PlutoUI
	using Plots
	using Random
	using Distributions
end

# ╔═╡ 10a9374c-8e1a-11eb-10b4-fb8272dd46d1
md"""
## SIR Model Simulation

#### Intro

The **SIR Model** has been developed to simulate an epidemic over time. The model consists of a system of 4 differential equations that express the rates of change of 4 variables over time. The 4 variables are:

1. **S** - the susceptibles of getting the infection
2. **I** - the infected
3. **R** - the recovered from the infection
4. **D** - the number of dead

#### The Model
Here follow the 4 equations that govern the model dynamics:

${ \begin{align*}
\frac{\mathrm{d}S}{\mathrm{d}t} & = -\beta \cdot \frac{S \cdot I}{N} + \delta \cdot R \\
\frac{\mathrm{d}I}{\mathrm{d}t} & = \beta \cdot \frac{S \cdot I}{N} + \gamma \cdot I \\
\frac{\mathrm{d}R}{\mathrm{d}t} & = \gamma \cdot (1-\rho) \cdot I - \delta \cdot R \\
\frac{\mathrm{d}D}{\mathrm{d}t} & = \gamma \cdot \rho \cdot I 
\end{align*} }$

#### Model simulation
The following is a simulation of the model described above. First choose the model parameters
"""

# ╔═╡ 43970fa0-8e1b-11eb-3052-1701966a4478
begin
	T = 38*30 	#period of 300 daysplot
	Δt = 1/400 	#time interval of 6 hours (1/4 of a day)
end;

# ╔═╡ 742d1fd0-907d-11eb-27d1-59ca2af59f97
parameter = Dict(
	"rate of infection" => 0.14,
	"rate of recovery" => 0.07,
	"rate of immuity loss" => 0.00,
	"rate of death of infected" => 0.001,
	"total population" => 6*10^2,
	"initial number of infected" => 10.0
);

# ╔═╡ d59335ac-913b-11eb-228f-273be8d7a86a
md"""
---
**Additional parameter and functions for simulations with two subgroups**
"""

# ╔═╡ 56b9df4c-913b-11eb-180a-eb784c2b5181
begin
	parameter_subpopulation = Dict(
		"rate of infection" => 0.14,
		"rate of recovery" => 0.035,
		"rate of immuity loss" => 0.0,
		"rate of death of infected" => 0.2,
		"total population" => 1*10^6,
		"initial number of infected" => 0.0
		)
	
	mixing_parameter = Dict(
		"rate of infection 1 -> 2" => 0.0,
		"rate of infection 2 -> 1" => 0.0
	)
end;

# ╔═╡ 6618e836-9080-11eb-39f4-21ac8e8af61a
md"""
---
**Additional Parameter and functions for simulations with an incidence dependent recovery rate with delay.**
"""

# ╔═╡ f8e23e02-9142-11eb-1700-6974da1a09c0
md"""
---
#### Worksheet 
Here all the calculations start. Set the initial state, choose the set of parameters and if needed the parameter changing function.
"""

# ╔═╡ e976a749-7fa0-4b81-b44d-461aec4527b8


# ╔═╡ da008842-9142-11eb-388c-d7e1f662d9db
md"""
---
#### Plots
"""

# ╔═╡ 93b1fcc4-9145-11eb-3d5e-4b9bd577a545
md"""
---
#### Model Functions

**Euler Algorithm to solve system. With constant or changing parameters** 
"""

# ╔═╡ 2680b820-8e32-11eb-28de-6f039360524b
begin
	function euler(f::Function,change_par!::Function,t0,tn,Δt,x₀,par)
		T = t0:Δt:tn
		F = Matrix{eltype(x₀)}(undef,(length(T),length(x₀)))
		F[1,:] = x₀
		for (n,t) ∈ enumerate(T[2:end])
			change_par!(par,t,F,n)
			F[n+1,:] = F[n,:] + f(t,F,n,par) . *1
		end
		return F
	end
	
	function no_change!(par,t,F,n) end
	
	euler(f::Function,t0,tn,Δt,x₀,par) = euler(f,no_change!,t0,tn,Δt,x₀,par)
end;

# ╔═╡ 58f42866-9146-11eb-29ea-87830bc8db66
md"""
**Defines the infinitesimal change of system**
"""

# ╔═╡ 3ab68bec-8e2f-11eb-203e-f9530bc39073
begin
	function SIR_step(t,x,β,δ,γ,ρ,N)
		S, I, R, D = x[1], x[2], x[3], x[4]
		nI = β*I*S/N*Δt
		return [
			-nI + δ*R*Δt,
			nI - γ*I,
			γ*(1-ρ)*I*Δt - δ*R*Δt,
			γ*ρ*I*Δt,
			nI
			]
	end
	SIR_step(t,x,par::Dict) = SIR_step(t,x,
				par["rate of infection"],par["rate of immuity loss"],
				par["rate of recovery"],par["rate of death of infected"],
				par["total population"]
			)
	SIR_step(t,H,n,par::Dict) = SIR_step(t,H[n,:],par)
end;

# ╔═╡ a7ac774a-913b-11eb-267b-8f9841c44cb4
begin
	function SIR_step_2D(t,H,n,par₁::Dict,par₂::Dict,par_mix::Dict)
		x₁,x₂ = view(H,n,1:4), view(H,n,5:8)	
		x = vcat(SIR_step(t,x₁,par₁),SIR_step(t,x₂,par₂))
		return x + SIR_interaction(t,x₁,x₂,
					par_mix["rate of infection 1 -> 2"],
					par_mix["rate of infection 2 -> 1"],
					par₁["total population"],par₂["total population"])
	end
	
	function SIR_interaction(t,x₁,x₂,β₃,β₄,N₁,N₂)
		S₁, I₁, S₂, I₂ = x₁[1], x₁[2], x₂[1], x₂[2]
		return [
			- β₄*S₁*I₂/N₁, β₄*S₁*I₂/N₁, 0, 0,
			- β₃*S₂*I₁/N₁, β₃*S₂*I₁/N₁, 0, 0
		]
	end
end;

# ╔═╡ de101aba-9145-11eb-2457-7b1831145eeb
md"""

**Supporting function for calculating the Incidence**
"""

# ╔═╡ 16a50518-906d-11eb-1e2b-3faa1765f30c
begin
	#Calculates at one time within the simulation
	function incidence(nI,t,pop_size;ndays=7,scale=100_000)
		#set starting time in range
		t₀ = max(t-round(Int,ndays/Δt),1)
		#return Incidence
		return (nI[t] - nI[t₀]) .* (scale/pop_size)
	end

	
	#Calculates for all times after simulation
	function incidence(nI,pop_size;ndays=7,scale=100_000)
		#rescale days
		ndays = round(Int,ndays/Δt)
		#sum new cases in first days
		I = nI[1:ndays]
		#iterate and add incidence
		for t₀ ∈ 1:length(nI)-ndays
			append!(I,nI[t₀+ndays] - nI[t₀])
		end
		#scale and return
		return I .* scale/pop_size
	end
end;

# ╔═╡ ce54bf76-906f-11eb-2480-75667d715789
begin
	#Additional parameter
	parameter_incidence_dependent_recovery_rate = merge(parameter,Dict(
			"high factor" => 2.8,
			"low incidence" => 40,
			"high incidence" => 200,
			"delay" => round(Int,5/Δt),
			"original recovery rate" => parameter["rate of recovery"]
			))
	
	#Dependency of rate of recovery on incidence
	function incidence_dependent_recovery_rate!(par,t,F,n)
		if n > par["delay"]
			i = incidence(view(F,:,5),n-par["delay"],par["total population"])
			if i > par["high incidence"]
				par["rate of recovery"] = par["original recovery rate"]*par["high factor"]
			end
			if i < par["low incidence"]
				par["rate of recovery"] = par["original recovery rate"]
			end
		end
	end
end;

# ╔═╡ c6bbb678-9080-11eb-2047-a3dcd34e8c77
begin
	N = parameter_incidence_dependent_recovery_rate["total population"]	
	I₀ = parameter_incidence_dependent_recovery_rate["initial number of infected"]
	
	x₀ = [
		N-I₀,
		I₀,
		0.0,
		0.0,
		0.0
	]
	
	result = euler(
		SIR_step,
		incidence_dependent_recovery_rate!,
		0.0,T,Δt,
		x₀,
		parameter_incidence_dependent_recovery_rate
		)
end;

# ╔═╡ 3d9185a6-908c-11eb-301c-358ed4a3eb37
let
	p = plot(framestlye=:zerolines,size=(680,450))

	plot!(p, 0.0:Δt:T,result[:,1:4],
		label = ["S" "I" "R" "D"],
		color = [:blue :red :green :yellow],
		xlabel = "Days",
		ylabel = "Number of Individuals"
	)
end

# ╔═╡ b2e6d55c-9142-11eb-3008-89e56b2af8db
let
	p = plot(framestlye=:zerolines,size=(680,450))

	plot!(p, 0.0:Δt:T,result[:,5],
		label = "I",
		color = :blue ,
		xlabel = "Days",
		ylabel = "Number of Cases"
	)
	
	p
end

# ╔═╡ d492c47e-908c-11eb-36bd-45fbb77d545e
let
	p = plot(framestlye=:zerolines,size=(680,450))

	I = incidence(result[:,5],parameter["total population"])

	plot!(p, 0.0:Δt:T,I,
		label = "I",
		color = :blue ,
		xlabel = "Days",
		ylabel = "7 day incidence per 100 000"
	)
	
	p
end

# ╔═╡ ceedca8e-9148-11eb-189c-a3a088520c2c
md"""
---
Package Management
"""

# ╔═╡ Cell order:
# ╟─10a9374c-8e1a-11eb-10b4-fb8272dd46d1
# ╠═43970fa0-8e1b-11eb-3052-1701966a4478
# ╠═742d1fd0-907d-11eb-27d1-59ca2af59f97
# ╟─d59335ac-913b-11eb-228f-273be8d7a86a
# ╠═56b9df4c-913b-11eb-180a-eb784c2b5181
# ╠═a7ac774a-913b-11eb-267b-8f9841c44cb4
# ╟─6618e836-9080-11eb-39f4-21ac8e8af61a
# ╠═ce54bf76-906f-11eb-2480-75667d715789
# ╟─f8e23e02-9142-11eb-1700-6974da1a09c0
# ╠═e976a749-7fa0-4b81-b44d-461aec4527b8
# ╠═c6bbb678-9080-11eb-2047-a3dcd34e8c77
# ╟─da008842-9142-11eb-388c-d7e1f662d9db
# ╠═3d9185a6-908c-11eb-301c-358ed4a3eb37
# ╟─d492c47e-908c-11eb-36bd-45fbb77d545e
# ╟─b2e6d55c-9142-11eb-3008-89e56b2af8db
# ╟─93b1fcc4-9145-11eb-3d5e-4b9bd577a545
# ╠═2680b820-8e32-11eb-28de-6f039360524b
# ╟─58f42866-9146-11eb-29ea-87830bc8db66
# ╠═3ab68bec-8e2f-11eb-203e-f9530bc39073
# ╟─de101aba-9145-11eb-2457-7b1831145eeb
# ╠═16a50518-906d-11eb-1e2b-3faa1765f30c
# ╟─ceedca8e-9148-11eb-189c-a3a088520c2c
# ╠═9f18ccd0-8e25-11eb-20e6-ab941f663550
