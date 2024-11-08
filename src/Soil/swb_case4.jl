function swb_case4(wa, IWS, pEc, pEs, s_tem, s_vod, soilpar, pftpar, fwet, zm, zgw)
  # INPUT:
  # wa      -- soil water content, 3 layers
  # IWS     -- total water enter into soil surface, mm
  # pEc     -- potential ET allocate to plant, mm
  # pEs     -- potential ET allocate to soil surface, mm
  # soilpar -- soil-related parameters
  # pftpar  -- plant-related parameters
  # wet     -- wetness indice
  # zm      -- soil layer depth, 3 layers
  # zgw     -- groundwater table depth, mm

  # Unsaturated depth in layer #1~3
  d1 = zm[1]
  d2 = zm[2]
  d3 = zm[3]

  wa1, wa2, wa3 = wa
  (; Ksat, θ_sat, θ_wp) = soilpar
  
  # ====== Water Supplement ====== #
  # Layer #1
  # Existed water column in the unsaturated zone #1
  wa1_unsat = wa1
  wc_s1 = d1 * wa1_unsat

  # Maximum water column in d1
  wc_m1 = d1 * θ_sat

  if wc_s1 + IWS >= wc_m1  # saturated
    wa1 = θ_sat  # current soil water content
    vw1 = wc_s1 + IWS - wc_m1  # exceeded water
  else  # unsaturated
    wa1 = wa1_unsat + IWS / d1  # soil water content in unsaturated zone
    vw1 = 0
  end

  # Layer #2
  # Existed water column in the unsaturated zone #2
  wa2_unsat = wa2
  wc_s2 = d2 * wa2_unsat

  # Maximum water column in d2
  wc_m2 = d2 * θ_sat

  if wc_s2 + vw1 >= wc_m2
    wa2 = θ_sat  # current soil water content
    vw2 = wc_s2 + vw1 - wc_m2  # exceeded water to layer #3
  else
    wa2 = wa2_unsat + vw1 / d2  # soil water content in unsaturated zone
    vw2 = 0  # no exceeded water
  end

  # Layer #3
  # Existed water column in the unsaturated zone #3
  wa3_unsat = wa3
  wc_s3 = d3 * wa3_unsat

  # Maximum water column in d3
  wc_m3 = d3 * θ_sat

  if wc_s3 + vw2 >= wc_m3
    wa3 = θ_sat  # current soil water content
    vw3 = wc_s3 + vw2 - wc_m3  # exceeded water
  else
    wa3 = wa3_unsat + vw2 / d3  # soil water content in unsaturated zone
    vw3 = 0  # no exceeded water
  end

  # ====== Water Consumption ====== #
  # Evapotranspiration #

  # Distributed the potential T to different layers
  Tr_p1, Tr_p2, Tr_p3 = pTr_partition(pEc, fwet, wa1, wa2, wa3, soilpar, pftpar, zm)

  # Calculate the moisture constrains for plant and soil in unsaturated zone
  # f_sm1, f_sm_s1 = swc_stress(wa1, pEc, soilpar, pftpar)
  # f_sm2, _ = swc_stress(wa2, pEc, soilpar, pftpar)
  # f_sm3, _ = swc_stress(wa3, pEc, soilpar, pftpar)

  f_sm1, f_sm_s1 = swc_stress(wa1, pEc, soilpar, pftpar)
  f_sm2, f_sm_s2 = swc_stress(wa2, pEc, soilpar, pftpar)
  f_sm3, f_sm_s3 = swc_stress(wa3, pEc, soilpar, pftpar)

  # f_sm1 = 1
  # f_sm_s1 = 1
  # f_sm2 = 1
  # f_sm_s2 = 1
  # f_sm3 = 1
  # f_sm_s3 = 1

  # Actual transpiration
  Tr1 = f_sm1 * s_vod * s_tem * Tr_p1
  Tr2 = f_sm2 * s_vod * s_tem * Tr_p2
  Tr3 = f_sm3 * s_vod * s_tem * Tr_p3
  Tr = Tr1 + Tr2 + Tr3

  # Actual soil evaporation
  Es = f_sm_s1 * pEs  # Only consider about the first layer

  # ====== Soil Water Drainage (Unsaturated Zone) ====== #
  # ---------------------------------------------------------------- Layer #1
  # Update the soil moisture after ET, layer #1
  Es = max(Es, 0)
  Tr1 = max(Tr1, 0)

  if wa1 > 0 && Es + Tr1 > d1 * wa1  # wilting point
    Tr1 = d1 * wa1 * Tr1 / (Tr1 + Es)
    Es = d1 * wa1 - Tr1
  end

  # Drainage from unsaturated zone, #1
  f1 = soil_drainage(wa1_unsat, θ_sat, Ksat, 0.048, 4.8)

  # Update the soil moisture after drainage, layer #1
  wa1 = (wa1 * d1 - f1 - Es - Tr1) / d1
  wa1 = clamp(wa1, 0, 1)

  # ---------------------------------------------------------------- Layer #2
  # Update the soil moisture after ET, layer #2
  Tr2 = clamp(Tr2, 0, d2 * (wa2 - θ_wp))

  # Drainage from unsaturated zone, #2
  f2 = soil_drainage(wa2_unsat, θ_sat, Ksat, 0.012, 1.2)

  # Update the soil moisture after drainage, layer #2
  wa2 = (wa2 * d2 + f1 - f2 - Tr2) / d2
  wa2 = clamp(wa2, 0, 1)  # > wwp 0

  if wa2 > θ_sat
    ff2 = max((wa2 - θ_sat) * d2, 0)  # extra water from upper layer
    wa2 = θ_sat
  else
    ff2 = 0
  end

  # ---------------------------------------------------------------- Layer #3
  # Update the soil moisture after ET, layer #3, unsat-zone
  Tr3 = clamp(Tr3, 0, d3 * (wa3 - θ_wp))

  # Drainage from unsaturated zone, #3
  f3 = soil_drainage(wa3_unsat, θ_sat, Ksat, 0.012, 1.2)

  # Update the soil moisture after drainage, layer #3
  wa3 = (wa3 * d3 + f2 + ff2 - f3 - Tr3) / d3

  wa3 = clamp(wa3, 0, 1)
  if wa3 > θ_sat
    ff3 = max((wa3 - θ_sat) * d3, 0)
    wa3_unsat = θ_sat
  else
    ff3 = 0
    wa3_unsat = wa3
  end

  # ====== The Groundwater Table Depth ====== #
  # Total water recharge to groundwater
  F1 = f3 + ff3 + vw3

  # Total transpiration from groundwater
  Tr_g = 0

  # R_sb groundwater discharge
  R_sb_max = 39  # mm day-1
  f = 1.25e-3  # mm-1
  R_sb = R_sb_max * exp(-f * zgw)

  # Variation of water stored in the saturated zone
  delta_w = F1 - Tr_g - R_sb

  # Changes in groundwater table depth
  delta_zgw = delta_w / 0.2  # specific yield as 0.2
  zgw = zgw - delta_zgw
  uex = 0  # excess water to soil surface, mm

  # Update soil moisture and groundwater table depth
  if zgw > zm[1] + zm[2] && zgw < zm[1] + zm[2] + zm[3]
    wa3 = (wa3_unsat * (zgw - zm[1] - zm[2]) + θ_sat * (zm[1] + zm[2] + zm[3] - zgw)) / zm[3]
  elseif zgw > zm[1] && zgw < zm[1] + zm[2]
    wa2 = (wa2 * (zgw - zm[1]) + θ_sat * (zm[1] + zm[2] - zgw)) / zm[2]
    wa3 = θ_sat
  elseif zgw > 0 && zgw < zm[1]
    wa1 = (wa1 * zgw + θ_sat * (zm[1] - zgw)) / zm[1]
    wa2 = θ_sat
    wa3 = θ_sat
  elseif zgw <= 0
    wa1 = θ_sat
    wa2 = θ_sat
    wa3 = θ_sat
    uex = -zgw * θ_sat  # excess water to soil surface, mm
  end

  # Updated soil water content
  wa = [wa1, wa2, wa3]
  zgw = max(0, zgw)
  return wa, zgw, Tr, Es, uex, Tr1, Tr2, Tr3, f_sm1, f_sm2, f_sm3, s_vod, s_tem,
  Tr_p1, Tr_p2, Tr_p3, f_sm_s1, f_sm_s2, f_sm_s3
end
