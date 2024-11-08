using SPAC, Ipaper, Dates
using RTableTools, DataFrames, NaNStatistics
using HydroTools
# using Ipaper, Ipaper.sf, ArchGDAL
# using GLMakie, MakieLayers

using Revise
activate

function init_param(; soil_type=3, PFTi=22)
  soilpar = get_soilpar(soil_type)
  pftpar = get_pftpar(PFTi)

  θ_sat = soilpar.θ_sat
  wa = ones(3) * θ_sat
  zgw = 0.0
  snowpack = 0.0
  state = State(; wa, zgw, snowpack)
  soilpar, pftpar, state
end

# begin
#   k = 1
#   x, y = st[k, [:lon, :lat]]
#   i, j = findnear(x, y, lon, lat)
#   soil_type = Soil[i, j]
#   soilpar = get_soilpar(soil_type)

#   topt = Float64(Topt[i, j])
#   soilpar, pftpar, state = init_param(;soil_type, PFTi=22)
# end
# Load necessary data
df = fread("data/dat_禹城_ERA5L_1979-2023.csv")

df_vpd = fread("E:/GitHub/SiTHv2_product/SiTHv2_kdd/data/EAR5L_VPD_01deg_1979-2023_yucheng.csv")
df_vpd.time = Date.(df_vpd.time .|> string, "yyyymmdd")
rename!(df_vpd, :time => :date)
rename!(df_vpd, :value => :VPD)

df_wind = fread("E:/GitHub/SiTHv2_product/SiTHv2_kdd/data/EAR5L_U10_01deg_1979-2023_yucheng.csv")
df_wind.time = Date.(df_wind.time .|> string, "yyyymmdd")
rename!(df_wind, :time => :date)
rename!(df_wind, :value => :Wind)

df = leftjoin(df, df_vpd, on = :date)
df = leftjoin(df, df_wind, on = :date)

dates = df.date

soilpar, pftpar, state = init_param()

inds = findall(year.(dates) .== 2010)
d = df[inds, :]
d.LAI = d.LAI |> drop_missing
d.VOD = d.VOD |> drop_missing
d.VPD = d.VPD |> drop_missing
d.Wind = d.Wind |> drop_missing
(; Rn, Pa, Prcp, Tavg, LAI, VOD, VPD, Wind) = d

Tas = deepcopy(Tavg) # Effective accumulated temperature
Tas[Tas.<0] .= 0 # Remove values less than 0
Tas = cumsum(Tas)

Gi = 0.4 .* Rn .* exp.(-0.5 .* LAI) # G_soil
s_VODi = (VOD ./ nanmaximum(VOD)) .^ 0.5 # VOD-stress

# s_VODi = (VOD ./ nanmaximum(VOD)) .^ 0.5 .* 0 .+ 1 # VOD-stress
topt = 24.0

ET, Tr, Es, Ei, Esb, SM, RF, GW, Tr1, Tr2, Tr3, case_num, f_sm1, f_sm2, f_sm3, s_vod, s_tem, Tr_p1, Tr_p2, Tr_p3, f_sm_s1, f_sm_s2, f_sm_s3, pEc, pEs, ET0 =
  SiTHv2_site(Rn, Tavg, Tas, Prcp, Pa, Gi, LAI, s_VODi, topt, soilpar, VPD, Wind, pftpar, state, false)

SM1 = SM[:, 1]
SM2 = SM[:, 2]
SM3 = SM[:, 3]

df_out = DataFrame(; ET, Tr, Es, Ei, Esb, SM1, SM2, SM3, RF, GW, Tr1, Tr2, Tr3, case_num, f_sm1, f_sm2, f_sm3, s_vod, s_tem, Tr_p1, Tr_p2, Tr_p3, f_sm_s1, f_sm_s2, f_sm_s3, pEc, pEs, ET0)
# fwrite(df_out, "data/Output_禹城_2010.csv")
fwrite(df_out, "data/Output_test_2010.csv")

# begin
#   using Plots
#   gr(framestyle = :box)
#   plot(ET)
# end
