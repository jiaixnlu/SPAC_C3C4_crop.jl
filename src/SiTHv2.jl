export SiTHv2_site

"""

# Model input
- Rn     : Net Radiation, W/m-2
- Ta     : Near surface air temperature, C
- Topt   : Optimal growth temperature for plant, C
- Pe     : Precipitation, mm day-1
- Pa     : Near surface air pressure, kPa
- s_VOD  : Constrains from VOD,[0,1]
- G      : Soil heat flux, W/m-2
- LAI    : Leaf area index
- soilpar: Parameters related to Soil types
- pftpar : Parameters related to PFTs
- wa     : Soil moisture (last step)
- zgw    : Groundwater table depth, mm
- snp    : Snow package (old), mm day-1

# Model output
- Et     : Total Evapotranspiration, mm day-1
- Tr     : Plant Transpiration, mm day-1
- Es     : Soil Evaporation, mm day-1
- Ei     : Intercepted Evaporation, mm day-1
- Esb    : Snow sublimation, mm day-1
- wa     : Soil moisture (three layers)
- srf    : Surface runoff, mm day-1
- zgw    : Groundwater table depth, mm
- snp    : Snow package (new), mm day-1

"""

function SiTHv2(Rn, Ta, Tas, Topt, P, Pa, s_VOD, G, LAI, soilpar, VPD, Wind, pftpar, state::State, i; K3=0.6, K4=0.6, K=0.6, Kc3=1.0, Kc4=1.0, Kc=1.0)
  (; wa, zgw, snowpack) = state
  pEc, pEs, ET0 = potentialET(Rn, G, LAI, Ta, Pa, VPD, Wind, i; K3=0.6, K4=0.6, K=0.6, Kc3=1.0, Kc4=1.0, Kc=1.0) # PET allocated to canopy and soil surface

  Ei, fwet, PE = interception(P, pEc, LAI, pftpar)  # Interception evaporation

  # Snow sublimation, snow melt
  state.snowpack, Esb, _, Pnet = snp_balance(PE, Ta, Tas, snowpack, pEs)

  srf, IWS, Vmax = runoff_up(Pnet, wa, zgw, ZM, soilpar)

  # Variables associated with soil water balance
  new_pEs = max(pEs - Esb, 0)
  Tr, Es, uex, Tr1, Tr2, Tr3, case_num, f_sm1, f_sm2, f_sm3, s_vod, s_tem, Tr_p1, Tr_p2, Tr_p3, f_sm_s1, f_sm_s2, f_sm_s3 = sw_balance(IWS, pEc, new_pEs, Ta, Topt, s_VOD, soilpar, pftpar, fwet, ZM, state)

  Et = Tr + Es + Ei + Esb # Total Evapotranspiration
  srf += uex
  return Et, Tr, Es, Ei, Esb, srf, Pnet, IWS, Vmax, Tr1, Tr2, Tr3, case_num, f_sm1, f_sm2, f_sm3, s_vod, s_tem, Tr_p1, Tr_p2, Tr_p3, f_sm_s1, f_sm_s2, f_sm_s3, pEc, pEs, ET0
end

function _run_model(Rn::T, Ta::T, Tas::T, Prcp::T, Pa::T, G::T, LAI::T, s_VOD::T,
  Top::FT, soilpar, VPD::T, Wind::T, pftpar, state::State, output; K3=0.6, K4=0.6, K=0.6, Kc3=1.0, Kc4=1.0, Kc=1.0) where {FT<:Real,T<:AbstractVector{FT}}
  ntime = size(Rn, 1)

  for i in 1:ntime
    
    # _Kc3 = 160 <= i <= 200 ? 0.7 : 1.0
    # _Kc4 = 160 <= i <= 200 ? 0.7 : 1.0
    # _Kc = 180 <= i <= 280 ? 1.0 : 1.0
    # Et, Tr, Es, Ei, Esb, srf, _, _, _, Tr1, Tr2, Tr3, case_num, f_sm1, f_sm2, f_sm3, s_vod, s_tem, Tr_p1, Tr_p2, Tr_p3, f_sm_s1, f_sm_s2, f_sm_s3, pEc, pEs, ET0 = SiTHv2(Rn[i], Ta[i], Tas[i], Top, Prcp[i], Pa[i], s_VOD[i], G[i], LAI[i], soilpar, VPD[i], Wind[i], pftpar, state; Kc=_Kc)

    Et, Tr, Es, Ei, Esb, srf, _, _, _, Tr1, Tr2, Tr3, case_num, f_sm1, f_sm2, f_sm3, s_vod, s_tem, Tr_p1, Tr_p2, Tr_p3, f_sm_s1, f_sm_s2, f_sm_s3, pEc, pEs, ET0 = SiTHv2(Rn[i], Ta[i], Tas[i], Top, Prcp[i], Pa[i], s_VOD[i], G[i], LAI[i], soilpar, VPD[i], Wind[i], pftpar, state, i; K3=0.6, K4=0.6, K=0.6, Kc3=1.0, Kc4=1.0, Kc=1.0)

    output[:ETs][i] = Et
    output[:Trs][i] = Tr
    output[:Ess][i] = Es
    output[:Eis][i] = Ei
    output[:Esbs][i] = Esb
    output[:SM][i, :] = state.wa
    output[:RF][i] = srf
    output[:GW][i] = state.zgw
    output[:Tr1][i] = Tr1
    output[:Tr2][i] = Tr2
    output[:Tr3][i] = Tr3
    output[:case_num][i] = case_num
    output[:f_sm1][i] = f_sm1
    output[:f_sm2][i] = f_sm2
    output[:f_sm3][i] = f_sm3
    output[:s_vod][i] = s_vod
    output[:s_tem][i] = s_tem
    output[:Tr_p1][i] = Tr_p1
    output[:Tr_p2][i] = Tr_p2
    output[:Tr_p3][i] = Tr_p3
    output[:f_sm_s1][i] = f_sm_s1
    output[:f_sm_s2][i] = f_sm_s2
    output[:f_sm_s3][i] = f_sm_s3
    output[:pEc][i] = pEc
    output[:pEs][i] = pEs
    output[:ET0][i] = ET0
  end
  return output
end


"""
## Arguments
- `Top`   : optimal growth temperature for plant, degC
- `state`: 
  + `wa`    : Soil moisture (last step)
  + `zgw`   : groundwater table depth, mm
  + `snp`   : Snow package (old), mm day-1
- `spinfg`: spin-up flag, 1 for spin-up, 0 for normal calculation. 循环重复跑100次。
"""
function SiTHv2_site(Rn, Ta, Tas, Prcp, Pa, G, LAI, s_VOD,
  Top, soilpar, VPD, Wind, pftpar, state, spin::Bool=false; K3=0.6, K4=0.6, K=0.6, Kc3=1.0, Kc4=1.0, Kc=1.0)

  ntime = length(Rn)
  output = Dict(
    :ETs => zeros(ntime),
    :Trs => zeros(ntime),
    :Ess => zeros(ntime),
    :Eis => zeros(ntime),
    :Esbs => zeros(ntime),
    :SM => zeros(ntime, 3),
    :RF => zeros(ntime),
    :GW => zeros(ntime),
    :Tr1 => zeros(ntime),
    :Tr2 => zeros(ntime),
    :Tr3 => zeros(ntime),
    :case_num => zeros(ntime),
    :f_sm1 => zeros(ntime),
    :f_sm2 => zeros(ntime),
    :f_sm3 => zeros(ntime),
    :s_vod => zeros(ntime),
    :s_tem => zeros(ntime),
    :Tr_p1 => zeros(ntime),
    :Tr_p2 => zeros(ntime),
    :Tr_p3 => zeros(ntime),
    :f_sm_s1 => zeros(ntime),
    :f_sm_s2 => zeros(ntime),
    :f_sm_s3 => zeros(ntime),
    :pEc => zeros(ntime),
    :pEs => zeros(ntime),
    :ET0 => zeros(ntime)
  )
  if spin == 1 # spin-up
    for k in 1:100 # set the spin-up time (100 years)
      output = _run_model(Rn, Ta, Tas, Prcp, Pa, G, LAI, s_VOD, Top, soilpar, VPD, Wind, pftpar, state, output; K3, K4, K, Kc3, Kc4, Kc)
    end
  else
    output = _run_model(Rn, Ta, Tas, Prcp, Pa, G, LAI, s_VOD, Top, soilpar, VPD, Wind, pftpar, state, output; K3, K4, K, Kc3, Kc4, Kc)
  end

  ETs = output[:ETs]
  Trs = output[:Trs]
  Ess = output[:Ess]
  Eis = output[:Eis]
  Esbs = output[:Esbs]
  SM = output[:SM]
  RF = output[:RF]
  GW = output[:GW]
  Tr1 = output[:Tr1]
  Tr2 = output[:Tr2]
  Tr3 = output[:Tr3]
  case_num = output[:case_num]
  f_sm1 = output[:f_sm1]
  f_sm2 = output[:f_sm2]
  f_sm3 = output[:f_sm3]
  s_vod = output[:s_vod]
  s_tem = output[:s_tem]
  Tr_p1 = output[:Tr_p1]
  Tr_p2 = output[:Tr_p2]
  Tr_p3 = output[:Tr_p3]
  f_sm_s1 = output[:f_sm_s1]
  f_sm_s2 = output[:f_sm_s2]
  f_sm_s3 = output[:f_sm_s3]
  pEc = output[:pEc]
  pEs = output[:pEs]
  ET0 = output[:ET0]
  return ETs, Trs, Ess, Eis, Esbs, SM, RF, GW, Tr1, Tr2, Tr3, case_num, f_sm1, f_sm2, f_sm3, s_vod, s_tem, Tr_p1, Tr_p2, Tr_p3, f_sm_s1, f_sm_s2, f_sm_s3, pEc, pEs, ET0
end
