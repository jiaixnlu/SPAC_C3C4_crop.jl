function swb_case0(wa, IWS, pEc, pEs, s_tem, s_vod, soilpar, pftpar, fwet, zm, zgw)
  # INPUT:
  # wa      -- soil water content, 3 layers
  # IWS     -- total water entering into soil surface (mm)
  # pEc     -- potential ET allocated to plants (mm)
  # pEs     -- potential ET allocated to soil surface (mm)
  # soilpar -- soil-related parameters
  # pftpar  -- plant-related parameters
  # wet     -- wetness index
  # zm      -- soil layer depth (3 layers)
  # zgw     -- groundwater table depth (mm)

  # Old soil water content in layer 1-3
  wa1, wa2, wa3 = wa
  (; θ_sat, θ_fc) = soilpar

  # ====== Water consumption ====== #
  # Evapotranspiration
  Tr_p1, Tr_p2, Tr_p3 = pTr_partition(pEc, fwet, wa1, wa2, wa3, soilpar, pftpar, zm)

  # Actual transpiration
  Tr1 = s_vod * s_tem * Tr_p1
  Tr2 = s_vod * s_tem * Tr_p2
  Tr3 = s_vod * s_tem * Tr_p3
  Tr = Tr1 + Tr2 + Tr3

  # Actual soil evaporation (only first layer)
  Es = pEs

  # Groundwater recharge and discharge calculations
  F1 = IWS
  Tr_g = Tr
  R_sb_max = 39  # maximum groundwater discharge (mm day-1)
  f = 1.25e-3    # discharge decay rate (mm-1)
  R_sb = R_sb_max * exp(-f * zgw)

  # Variation of water stored in the saturated zone
  delta_w = F1 - Tr_g - R_sb

  # Change in groundwater table depth
  delta_zgw = delta_w / (θ_sat - θ_fc)
  zgw -= delta_zgw
  uex = 0  # excess water to soil surface

  # Update soil moisture and groundwater table depth
  if zgw > zm[1] + zm[2] + zm[3]
    wa1 = θ_fc
    wa2 = θ_fc
    wa3 = θ_fc
  elseif zgw > zm[1] + zm[2] && zgw <= zm[1] + zm[2] + zm[3]
    wa1 = θ_fc
    wa2 = θ_fc
    wa3 = (θ_fc * (zgw - zm[1] - zm[2]) + θ_sat * (zm[1] + zm[2] + zm[3] - zgw)) / zm[3]
  elseif zgw > zm[1] && zgw <= zm[1] + zm[2]
    wa1 = θ_fc
    wa2 = (θ_fc * (zgw - zm[1]) + θ_sat * (zm[1] + zm[2] - zgw)) / zm[2]
    wa3 = θ_sat
  elseif zgw > 0 && zgw <= zm[1]
    wa1 = (θ_fc * zgw + θ_sat * (zm[1] - zgw)) / zm[1]
    wa2 = θ_sat
    wa3 = θ_sat
  elseif zgw <= 0
    wa1 = θ_sat
    wa2 = θ_sat
    wa3 = θ_sat
    uex = -zgw * θ_fc  # excess water to soil surface
  end

  # Updated soil water content
  wa = [wa1, wa2, wa3]
  zgw = max(0, zgw)  # ensure groundwater table depth is non-negative

  f_sm1 = 1
  f_sm2 = 1
  f_sm3 = 1
  f_sm_s1 = 1
  f_sm_s2 = 1 
  f_sm_s3 = 1
  return wa, zgw, Tr, Es, uex, Tr1, Tr2, Tr3, f_sm1, f_sm2, f_sm3, s_vod, s_tem, Tr_p1, Tr_p2, Tr_p3, f_sm_s1, f_sm_s2, f_sm_s3
end
