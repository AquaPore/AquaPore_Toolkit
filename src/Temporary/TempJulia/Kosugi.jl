	import SpecialFunctions: erfc, erfcinv


function RUN_KOSUGI()
	function KG(Ψ; θs, θr, σ)
		Ψm =  10. * exp(σ * 3.0)
		θ_Mat = 0.5 * (θs - θr) * erfc((log(Ψ / Ψm)) / (σ * √2.0)) + θr
	return θ_Mat
	end


	Ψ₀ = [0.0, 500.0, 1000.0, 2000.0,  4000.0,  10000.0,  150000.0]
	N = length(Ψ₀)
	θ =fill(0.0,N)

	Rf=0.0

	for iZ=1:N
		θ[iZ] = KG(Ψ₀[iZ], θs=0.1/(1-Rf), θr=0.0/(1-Rf), σ=3.)
	end 

	println(Ψ₀)
	println(θ)

return nothing
end

RUN_KOSUGI()





