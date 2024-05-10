



function  θΨ_MacMat(Ψ, θs, θsMacMat, θr, Ψm, σ, ΨmMac, ΨmacMat, σMac)
		if Ψ₁ ≤ ΨmacMat
			return θ_Mac =  max(θs - θsMacMat, 0.0) * 0.5 * (erfc((log(Ψ / ΨmMac)) / (σMac * √2.0)) - (min(Ψ / ΨmacMat, 1.0)) * erfc((log(ΨmacMat / ΨmMac)) / (σMac * √2.0))) + θsMacMat
		else
			return θ_Mat = 0.5 * (θsMacMat - θr) * erfc((log( max(Ψ₁ - ΨmacMat, eps(100.0)) / Ψm)) / (σ * √2.0)) + θr
		end
end 
