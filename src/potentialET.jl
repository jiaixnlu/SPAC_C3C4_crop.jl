const Cp = 1013  # Specific heat (J kg-1 C-1)
const ϵ = 0.622  # e (unitless) is the ratio of molecular weight of water to dry air (equal to 0.622)

"""
Potential ET partition

## INPUT:
- Ta    : air temperature, C
- Rn    : average daily net radiation, W/m^2
- Pa    : atmospheric pressure, kPa
- LAI   : Leaf area index, 1
- G     : Soil heat flux, W/m^2

## OUTPUT
- pEc   : potential Transpiration, mm/day
- pEs   : potential Soil evaporation, mm/day
"""
function potentialET(Rn::T, G::T, LAI::T, Ta::T, Pa::T, VPD::T, Wind::T, i::Int; K3::T=0.6, K4::T=0.6, K::T=0.6, Kc3::T=1.0, Kc4::T=1.0, Kc::T=1.0) where {T<:Real}

  Kc3 = 0.75
  Kc4 = 0.9

  if 32 <= i <= 166
    Kc = Kc3
  elseif 196 < i <= 288
    Kc = Kc4
  else
    Kc = Kc
  end

  K3 = 0.6
  K4 = 0.6
  K = 0.6
  if 32 <= i <= 166
    k = K3
  elseif 196 < i <= 288
    k = K4
  else
    k = K
  end

  # k = 0.6  # the empirical extinction coefficient set as 0.6
  # Radiation located into soil and canopy, separately
  Rns = exp(-k * LAI) * Rn
  Rnc = Rn - Rns

  λ = cal_lambda(Ta)                   # [J/kg]
  Δ = cal_slope(Ta)                    # [kPa/degC]
  γ = (Cp * Pa) / (ϵ * λ)              # Psychrometric constant (kPa/degC)

  Eeq_canopy = Δ / (Δ + γ) * Rnc
  Eeq_canopy = Eeq_canopy * 86400 / λ

  # Eeq_soil = Δ / (Δ + γ) * (Rns - G)
  # Eeq_soil = Eeq_soil * 86400 / λ

  Eeq_total = Δ / (Δ + γ) * Rn
  Eeq_total = Eeq_total * 86400 / λ

  # Potential Transpiration and Soil evaporation, mm/day
  # pEc = ET0_PT1972(Rnc, Δ, γ, λ) * Kc
  # pEc = ET0_PT1972(Rnc, Δ, γ, λ) * Kc
  pEc = ET0_Penman48(Δ, γ, λ, Eeq_canopy, VPD, Wind) * Kc
  pEs = ET0_PT1972(Rns - G, Δ, γ, λ)
  # pEs = ET0_Penman48(Δ, γ, λ, Eeq_soil, VPD, Wind)
  # ET0 = ET0_PT1972(Rn, Δ, γ, λ)
  ET0 = ET0_Penman48(Δ, γ, λ, Eeq_total, VPD, Wind)
  return pEc, pEs, ET0
end


# Rn: W m-2
function ET0_PT1972(Rn::T, Δ::T, γ::T, λ::T) where {T<:Real}
  α = 1.26  # PT coefficient for water saturated surface
  ET0 = α * Rn * Δ / (Δ + γ)
  ET0 = max(ET0, 0)
  W2mm(ET0, λ)
end

# λ: latent heat of vaporization (J kg-1)
W2mm(w::T, λ::T) where {T<:Real} = w * 86400 / λ

# λ: latent heat of vaporization (J kg-1)
cal_lambda(Ta::T) where {T<:Real} = 2.501e6 - 2361 * Ta

# cal_gamma(Pa, λ) = (Cp * Pa) / (ϵ * λ)
cal_es(Ta::T) where {T<:Real} = 0.6108 * exp((17.27 * Ta) / (Ta + 237.3))

cal_slope(Ta::T) where {T<:Real} = 4098 * cal_es(Ta) / ((Ta + 237.3)^2)


function ET0_Penman48(Δ::T, γ::T, λ::T, Eeq::T, VPD::T, wind::T; z_wind=2) where {T<:Real}

  U2= cal_U2(wind, z_wind)
  Evp = γ / (Δ + γ) * 6.43 * (1 + 0.536 * U2) * VPD / λ * 1e6
  ET0 = Evp + Eeq
  ET0
end

function cal_U2(wind, z_wind = 10)
  if z_wind == 2
    return wind
  end
  return wind * 4.87 / log(67.8 * z_wind - 5.42)
end

# # 假设输入值
# Rn = 200.0      # W/m²
# Ta = 25.0       # °C
# Pa = 101.3      # kPa
# VPD = 1.5       # kPa
# wind = 3.0      # m/s


# λ = cal_lambda(Ta)                   # [J/kg]  
# Δ = cal_slope(Ta)                    # [kPa/degC]
# γ = (Cp * Pa) / (ϵ * λ)              # Psychrometric constant (kPa/degC)

# ET0_PT1972(Rn, Δ, γ, λ)

# Eeq = Δ / (Δ + γ) * Rn
# Eeq = Eeq * 86400 / λ

# U2 = cal_U2(wind, z_wind)

# Evp = γ / (Δ + γ) * 6.43 * (1 + 0.536 * U2) * VPD / λ * 1e6


