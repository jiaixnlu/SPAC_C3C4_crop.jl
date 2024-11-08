function swb_case3(wa, IWS, pEc, pEs, s_tem, s_vod, soilpar, pftpar, fwet, zm, zgw)
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
  d3 = zgw - zm[1] - zm[2]
  
  wa1, wa2, wa3 = wa
  (; Ksat, θ_sat, θ_fc, θ_wp) = soilpar
  
  # ====== Water Supplement ====== #
  # Layer #1
  wa1_unsat = wa1
  wc_s1 = d1 * wa1_unsat
  wc_m1 = d1 * θ_sat

  if wc_s1 + IWS >= wc_m1
    wa1 = θ_sat
    vw1 = wc_s1 + IWS - wc_m1
  else
    wa1 = wa1_unsat + IWS / d1
    vw1 = 0
  end

  # Layer #2
  wa2_unsat = wa2
  wc_s2 = d2 * wa2_unsat
  wc_m2 = d2 * θ_sat

  if wc_s2 + vw1 >= wc_m2
    wa2 = θ_sat
    vw2 = wc_s2 + vw1 - wc_m2
  else
    wa2 = wa2_unsat + vw1 / d2
    vw2 = 0
  end

  # Layer #3
  wa3_unsat = (wa3 * zm[3] - θ_sat * (zm[3] - d3)) / d3
  wc_s3 = d3 * wa3_unsat
  wc_m3 = d3 * θ_sat

  if wc_s3 + vw2 >= wc_m3
    wa3 = θ_sat
    vw3 = wc_s3 + vw2 - wc_m3
  else
    wa3_unsat = wa3_unsat + vw2 / d3
    wa3 = (wa3_unsat * d3 + θ_sat * (zm[3] - d3)) / zm[3]
    vw3 = 0
  end

  # ====== Water Consumption ====== #
  Tr_p1, Tr_p2, Tr_p3 = pTr_partition(pEc, fwet, wa1, wa2, wa3, soilpar, pftpar, zm)
  Tr_p3_u = Tr_p3 * (d3 * wa3_unsat) / (d3 * wa3_unsat + (zm[3] - d3) * θ_sat)
  Tr_p3_g = Tr_p3 * ((zm[3] - d3) * θ_sat) / (d3 * wa3_unsat + (zm[3] - d3) * θ_sat)

  # f_sm1, f_sm_s1 = swc_stress(wa1, pEc, soilpar, pftpar)
  # f_sm2, _ = swc_stress(wa2, pEc, soilpar, pftpar)
  # f_sm3, _ = swc_stress(wa3_unsat, pEc, soilpar, pftpar)

  f_sm1, f_sm_s1 = swc_stress(wa1, pEc, soilpar, pftpar)
  f_sm2, f_sm_s2 = swc_stress(wa2, pEc, soilpar, pftpar)
  f_sm3, f_sm_s3 = swc_stress(wa3_unsat, pEc, soilpar, pftpar)

  # f_sm1 = 1
  # f_sm_s1 = 1
  # f_sm2 = 1
  # f_sm_s2 = 1
  # f_sm3 = 1
  # f_sm_s3 = 1

  Tr1 = f_sm1 * s_vod * s_tem * Tr_p1
  Tr2 = f_sm2 * s_vod * s_tem * Tr_p2
  Tr3_u = f_sm3 * s_vod * s_tem * Tr_p3_u
  Tr3_g = s_vod * s_tem * Tr_p3_g
  Tr3 = Tr3_u + Tr3_g
  Tr = Tr1 + Tr2 + Tr3

  Es = f_sm_s1 * pEs

  # ====== Soil Water Drainage (Unsaturated Zone) ====== #
  Es = max(Es, 0)
  Tr1 = max(Tr1, 0)

  if wa1 > 0 && Es + Tr1 > d1 * wa1
    Tr1 = d1 * wa1 * Tr1 / (Tr1 + Es)
    Es = d1 * wa1 - Tr1
  end

  f1 = soil_drainage(wa1_unsat, θ_sat, Ksat, 0.048, 4.8)
  wa1 = (wa1 * d1 - f1 - Es - Tr1) / d1
  wa1 = max(wa1, 0)

  Tr2 = clamp(Tr2, 0, d2 * (wa2 - θ_wp))
  f2 = soil_drainage(wa2_unsat, θ_sat, Ksat, 0.012, 1.2)
  wa2 = (wa2 * d2 + f1 - f2 - Tr2) / d2
  wa2 = clamp(wa2, 0, 1)

  if wa2 > θ_sat
    ff2 = max((wa2 - θ_sat) * d2, 0)
    wa2 = θ_sat
  else
    ff2 = 0
  end

  Tr3_u = max(Tr3_u, 0)
  Tr3_u = min(Tr3_u, d3 * (wa3_unsat - θ_wp))
  f3 = soil_drainage(wa3_unsat, θ_sat, Ksat, 0.012, 1.2)
  wa3_unsat = (wa3_unsat * d3 + f2 + ff2 - f3 - Tr3_u) / d3
  wa3_unsat = max(wa3_unsat, 0)
  wa3_unsat = min(wa3_unsat, 1)

  if wa3_unsat > θ_sat
    ff3 = max((wa3_unsat - θ_sat) * d3, 0)
    wa3_unsat = θ_sat
  else
    ff3 = 0
  end

  # ====== The Groundwater Table Depth ====== #
  F1 = f3 + ff3 + vw3
  Tr_g = Tr3_g
  R_sb_max = 39
  f = 1.25e-3
  R_sb = R_sb_max * exp(-f * zgw)
  delta_w = F1 - Tr_g - R_sb
  delta_zgw = delta_w / (θ_sat - (wa1 + wa2 + wa3_unsat) / 3)
  zgw = zgw - delta_zgw
  uex = 0

  if zgw > zm[1] + zm[2] + zm[3]
    wa3 = (wa3_unsat * d3 + θ_fc * (zm[3] - d3)) / zm[3]
  elseif zgw > zm[1] + zm[2] && zgw < zm[1] + zm[2] + zm[3]
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
    uex = -zgw * θ_sat
  end

  wa = [wa1, wa2, wa3]
  zgw = max(0, zgw)
  return wa, zgw, Tr, Es, uex, Tr1, Tr2, Tr3, f_sm1, f_sm2, f_sm3, s_vod, s_tem, Tr_p1, Tr_p2, Tr_p3, f_sm_s1, f_sm_s2, f_sm_s3
end
