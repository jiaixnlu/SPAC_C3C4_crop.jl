"""
Snowpack balance

- snowpack : available snow storage
- snowmelt : snow melt
- Esb      : snow sublimation
"""
function snp_balance(Prcp, Ta, Tas, snowpack, pEs)
  # Esnow_emp = 0.84 * (0.864 * (7.093 * Ta + 28.26)) / (Ta^2 - 3.593 * Ta + 5.175)
  Esnow = pEs  # Simple equivalent (Needs further development)

  # only snowfall occurs at Ta below zero
  if Ta <= 0
    # Add new snowfall, Ta <= 0
    newsnow = Prcp
    snowpack += newsnow

    snowmelt = 0  # snow melt
    Esb = clamp(snowpack, 0, Esnow)  # real snow sublimation

    Pnet = 0  # net Precipitation into soil surface
    snowpack -= Esb  # new snowpack
  else
    # real snow sublimation
    Esb = clamp(snowpack, 0, Esnow)

    snowpack -= Esb
    # snow melt, Ta > 0
    snowmelt_x = (1.5 + 0.007 * Prcp) * Tas  # Tas, accumulated Ta > 0
    snowmelt = min(snowpack, snowmelt_x)
    snowpack -= snowmelt

    Pnet = max(0, Prcp + snowmelt)  # net water into soil surface
  end
  return snowpack, Esb, snowmelt, Pnet
end
