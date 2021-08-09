### A Pluto.jl notebook ###
# v0.15.1

using Markdown
using InteractiveUtils

# ╔═╡ 9f18ccd0-8e25-11eb-20e6-ab941f663550
begin
	import Pkg
    Pkg.activate(".")
    Pkg.add("PlutoUI")
	Pkg.add("Plots")
	Pkg.add("Distributions")
	
	using PlutoUI
	using Plots
	using Distributions
end

# ╔═╡ 10a9374c-8e1a-11eb-10b4-fb8272dd46d1
md"""
## Random SIRD Model with two Populations
#### Intro
The **SIRD Model** has been developed to simulate an epidemic over time. The model consists of a system of 4 differential equations that express the rates of change of 4 variables over time. The 4 variables are:
${ \begin{align*}
	S & - \text{the susceptibles of getting the infection} \\
	I & - \text{the infected} \\
	R & - \text{the recovered from the infection} \\
	D & - \text{the number of dead}
\end{align*} }$
#### The Model
Define the parameters of the dynamics as follows:
${ \begin{align*}
	\beta & \to \text{ rate of infection} \\
	\gamma & \to \text{ rate of recovery} \\
	\delta & \to \text{ rate of immunity loss} \\
	\rho & \to \text{ mortality}
\end{align*} }$
Then the following 4 equations govern the dynamis of the model:
${ \begin{align*}
\frac{\mathrm{d}S}{\mathrm{d}t} & = -\beta \cdot \frac{S \cdot I}{N} + \delta \cdot R \\
\frac{\mathrm{d}I}{\mathrm{d}t} & = \beta \cdot \frac{S \cdot I}{N} - (\gamma+\rho) \cdot I \\
\frac{\mathrm{d}R}{\mathrm{d}t} & = \gamma \cdot I - \delta \cdot R \\
\frac{\mathrm{d}D}{\mathrm{d}t} & = \rho \cdot I 
\end{align*} }$
To model two mixing populations with independent parameters double the system and add resp. substract the mixing term 
${ \pm \beta_{ij} \frac{S_j\cdot I_i}{N_i + N_j} }$
to the number of infected resp. susceptible individuals. This term models the infection from Population \\( i \\) to Population \\( j \\), for \\( i,j = 1,2\\).
Finally to randomize the dynamics of the system in every step, all quantities that move between within the system get drawn from a Poisson distribution with the respective rate as mean.
#### Model simulation
The following is a simulation of the model described above. First choose the model parameters
"""

# ╔═╡ 43970fa0-8e1b-11eb-3052-1701966a4478
begin
	T = 3000 	#period of 300 daysplot
	Δt = 1/40 	#time interval of 6 hours (1/4 of a day)
end;

# ╔═╡ 742d1fd0-907d-11eb-27d1-59ca2af59f97
parameter = Dict(
	"rate of infection" => 0.14,
	"rate of recovery" => 0.07,
	"rate of immunity loss" => 0.01,
	"rate of death of infected" => 0.001,
	"total population" => 1*10^8,
	"initial number of infected" => 10
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
		"rate of immunity loss" => 0.02,
		"rate of death of infected" => 0.01,
		"total population" => 1*10^2,
		"initial number of infected" => 0.0
		)
	
	mixing_parameter = Dict(
		"rate of infection 1 -> 2" => 0.00014,
		"rate of infection 2 -> 1" => 0.14
	)
end;

# ╔═╡ f8e23e02-9142-11eb-1700-6974da1a09c0
md"""
---
#### Worksheet 
Here all the calculations start. Set the initial state, choose the set of parameters and if needed the parameter changing function.
"""

# ╔═╡ da008842-9142-11eb-388c-d7e1f662d9db
md"""
---
#### Plots
"""

# ╔═╡ deb93731-e412-47a6-b615-3ca7d82052dd


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
			F[n+1,:] = max.(F[n,:] + f(t,F,n,par) .*1,0)
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
		nI = rand(Poisson(β*I*S/N*Δt))
		nIL = rand(Poisson(δ*R*Δt))
		nR = rand(Poisson(γ*(1-ρ)*I*Δt))
		nD = rand(Poisson(ρ*γ*I*Δt))
		return [
			-nI + nIL,
			nI - nR - nD,
			nR - nIL,
			nD
			]
	end
	
	SIR_step(t,x,par::Dict) = SIR_step(t,x,
				par["rate of infection"],par["rate of immunity loss"],
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
					par₁["total population"]+par₂["total population"])
	end
	
	SIR_step_2D(t,H,n,par) = SIR_step_2D(t,H,n,par[1],par[2],par[3])
	
	function SIR_interaction(t,x₁,x₂,β₃,β₄,N)
		S₁, I₁, S₂, I₂ = x₁[1], x₁[2], x₂[1], x₂[2]
		nI₁ = rand(Poisson(β₄*S₁*I₂/N*Δt))
		nI₂ = rand(Poisson(β₃*S₂*I₁/N*Δt))
		return [
			-nI₁, nI₁, 0, 0,
			-nI₂, nI₂, 0, 0
		]
	end
end;

# ╔═╡ c6bbb678-9080-11eb-2047-a3dcd34e8c77
begin
	N₁,N₂ = parameter["total population"], parameter_subpopulation["total population"]
	I₁,I₂ = parameter["initial number of infected"], parameter_subpopulation["initial number of infected"]
	
	x₀ = [
		N₁-I₁,
		I₁,
		0.0,
		0.0,
		N₂-I₂,
		I₂,
		0.0,
		0.0
	]
	
	result = euler(
		SIR_step_2D,
		0.0,T,Δt,
		x₀,
		(parameter,parameter_subpopulation,mixing_parameter)
		)
end;

# ╔═╡ 3d9185a6-908c-11eb-301c-358ed4a3eb37
let
	p1 = plot(framestlye=:zerolines,title="Main Population")

	plot!(p1, 0.0:Δt:T,result[:,1:4],
		label = ["S" "I" "R" "D"],
		color = [:blue :red :green :yellow],
		xlabel = "Days",
		ylabel = "Number of Individuals"
	)
	
	p2 = plot(framestlye=:zerolines,title="Subpopulation")

	plot!(p2, 0.0:Δt:T,result[:,5:8],
		label = ["S" "I" "R" "D"],
		color = [:blue :red :green :yellow],
		xlabel = "Days",
		ylabel = "Number of Individuals"
	)
	
	plot(p1,p2,layout=(2,1),size=(680,680))
end

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
# ╟─f8e23e02-9142-11eb-1700-6974da1a09c0
# ╠═c6bbb678-9080-11eb-2047-a3dcd34e8c77
# ╟─da008842-9142-11eb-388c-d7e1f662d9db
# ╠═deb93731-e412-47a6-b615-3ca7d82052dd
# ╟─3d9185a6-908c-11eb-301c-358ed4a3eb37
# ╟─93b1fcc4-9145-11eb-3d5e-4b9bd577a545
# ╠═2680b820-8e32-11eb-28de-6f039360524b
# ╟─58f42866-9146-11eb-29ea-87830bc8db66
# ╠═3ab68bec-8e2f-11eb-203e-f9530bc39073
# ╠═a7ac774a-913b-11eb-267b-8f9841c44cb4
# ╟─de101aba-9145-11eb-2457-7b1831145eeb
# ╠═16a50518-906d-11eb-1e2b-3faa1765f30c
# ╟─ceedca8e-9148-11eb-189c-a3a088520c2c
# ╠═9f18ccd0-8e25-11eb-20e6-ab941f663550
