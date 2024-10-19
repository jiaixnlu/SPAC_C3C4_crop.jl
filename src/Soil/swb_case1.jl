function swb_case1(wa, IWS, pEc, pEs, s_tem, s_vod, soilpar, pftpar, wet, zm, zgw)
  # INPUT:
  # wa      -- soil water content, 3 layers
  # IWS     -- total water entering into the soil surface (mm)
  # pEc     -- potential ET allocated to plants (mm)
  # pEs     -- potential ET allocated to soil surface (mm)
  # soilpar -- soil-related parameters
  # pftpar  -- plant-related parameters
  # wet     -- wetness index
  # zm      -- soil layer depths, 3 layers
  # zgw     -- groundwater table depth (mm)

  # Unsaturated depth in layer #1
  d1 = zgw

  wa1, wa2, wa3 = wa

  ks = soilpar[1]  # hydraulic conductivity
  theta_sat = soilpar[3]  # saturated soil water content
  theta_fc = soilpar[5]  # field capacity
  wwp = soilpar[7]  # wilting point

  # ====== Water supplement ====== #
  # Layer #1 - Unsaturated zone
  wa1_unsat = (wa1 * zm[1] - theta_sat * (zm[1] - d1)) / d1  # Unsaturated region soil moisture
  wc_s1 = d1 * wa1_unsat  # Current water column (mm)
  wc_m1 = d1 * theta_sat  # Maximum water column (mm)

  if wc_s1 + IWS >= wc_m1
    wa1_unsat = theta_sat
    vw1 = wc_s1 + IWS - wc_m1  # Excess water
  else
    wa1_unsat += IWS / d1
    wa1 = (wa1_unsat * d1 + theta_sat * (zm[1] - d1)) / zm[1]
    vw1 = 0
  end

  # Layer #2 and #3 - Fully saturated

  # ====== Water consumption ====== #
  # Evapotranspiration
  Tr_p1, Tr_p2, Tr_p3 = pTr_partition(pEc, wa1, wa2, wa3, soilpar, pftpar, wet, zm)

  # Transpiration from unsaturated and saturated zones in layer #1
  Tr_p1_u = Tr_p1 * (d1 * wa1_unsat) / (d1 * wa1_unsat + (zm[1] - d1) * theta_sat)
  Tr_p1_g = Tr_p1 * ((zm[1] - d1) * theta_sat) / (d1 * wa1_unsat + (zm[1] - d1) * theta_sat)

  # Moisture constraints
  f_sm1, f_sm_s1 = swc_stress(wa1, soilpar, pEc, pftpar)

  # Actual transpiration
  Tr1_u = clamp(f_sm1 * s_vod * s_tem * Tr_p1_u, 0, d1 * (wa1_unsat - wwp))
  Tr1_g = s_vod * s_tem * Tr_p1_g
  Tr1 = Tr1_u + Tr1_g
  Tr2 = s_vod * s_tem * Tr_p2
  Tr3 = s_vod * s_tem * Tr_p3
  Tr = Tr1 + Tr2 + Tr3

  # Actual soil evaporation
  Es = f_sm_s1 * pEs
  Es_u = clamp(Es * (d1 * wa1_unsat) / (d1 * wa1_unsat + (zm[1] - d1) * theta_sat), 0, d1 * wa1_unsat)

  # ====== Soil water drainage (unsaturated zone) ====== #
  f1 = soil_drainage(wa1_unsat, theta_sat, ks, 0.048, 4.8)
  wa1_unsat = max((wa1_unsat * d1 - f1 - Es_u - Tr1_u) / d1, 0)

  # ====== Groundwater table depth update ====== #
  F1 = f1 + vw1  # Total water recharge to groundwater
  Tr_g = Tr1_g + Tr2 + Tr3  # Total transpiration from groundwater

  # Groundwater discharge
  R_sb_max = 39  # mm day-1
  f = 1.25e-3  # mm-1
  R_sb = R_sb_max * exp(-f * zgw)

  # Variation of water stored in the saturated zone
  delta_w = F1 - Tr_g - R_sb

  # Change in groundwater table depth
  delta_zgw = delta_w / (theta_sat - wa1_unsat)
  zgw -= delta_zgw
  uex = 0  # Excess water to the soil surface

  # Update soil moisture and groundwater table depth
  if zgw > zm[1] + zm[2] + zm[3]
    wa1 = (wa1_unsat * d1 + theta_fc * (zm[1] - d1)) / zm[1]
    wa2 = theta_fc
    wa3 = theta_fc
  elseif zgw > zm[1] + zm[2] && zgw <= zm[1] + zm[2] + zm[3]
    wa1 = (wa1_unsat * d1 + theta_fc * (zm[1] - d1)) / zm[1]
    wa2 = theta_fc
    wa3 = (theta_fc * (zgw - zm[1] - zm[2]) + theta_sat * (zm[1] + zm[2] + zm[3] - zgw)) / zm[3]
  elseif zgw > zm[1] && zgw <= zm[1] + zm[2]
    wa1 = (wa1_unsat * d1 + theta_fc * (zm[1] - d1)) / zm[1]
    wa2 = (theta_fc * (zgw - zm[1]) + theta_sat * (zm[1] + zm[2] - zgw)) / zm[2]
    wa3 = theta_sat
  elseif zgw > 0 && zgw <= zm[1]
    wa1 = (wa1_unsat * zgw + theta_sat * (zm[1] - zgw)) / zm[1]
    wa2 = theta_sat
    wa3 = theta_sat
  elseif zgw <= 0
    wa1 = theta_sat
    wa2 = theta_sat
    wa3 = theta_sat
    uex = -zgw * theta_fc  # Excess water to soil surface
  end

  # Updated soil water content
  wa = [wa1, wa2, wa3]
  zgw = max(0, zgw)

  return wa, zgw, Tr, Es, uex
end
