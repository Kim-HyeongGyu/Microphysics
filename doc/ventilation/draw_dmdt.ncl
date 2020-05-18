data1 = asciiread("./ON_dm_dt_nt1_nz1.txt", -1,"float")
data2 = asciiread("./OFF_dm_dt_nt1_nz1.txt",-1,"float")
data3 = asciiread("./r.txt",-1,"float")

x = ispan(1,100,1)

data = new((/3,dimsizes(x)/),"float")
data(0,:) = data1       ; ventilation ON
data(1,:) = data2       ; ventilation OFF
data(2,:) = data3       ; r 

Ttext = "T =  293.15 K"
Ptext = "P = 1011.94 hPa"

wks = gsn_open_wks("png","ventilation")

res = True
res@gsnFrame    = False
res@gsnMaximize = True
res@vpHeightF  = 0.5
res@vpWidthF   = 0.5

res@xyLineColors        = (/"blue","black"/)
res@xyDashPatterns      = (/0.,0./)
res@xyLineThicknesses   = (/3.,3./)

res@tiMainString    = "Ventilation Effect"
res@tiXAxisString   = "r [m]"
res@tiYAxisString   = "dm/dt [kg/s]"

res@trYMinF = 0.
res@trYMaxF = 1.5e-13

res@tmXTOn = False
res@tmYROn = False

res@xyXStyle    = "Log"

res@pmLegendDisplayMode = "Always"
res@pmLegendSide        = "Top"
res@pmLegendHeightF     = 0.04
res@pmLegendWidthF      = 0.15
res@lgPerimOn           = False
res@lgLabelFontHeightF  = 0.01
res@xyExplicitLegendLabels  = (/"   on","   off"/)
res@lgItemOrder         = (/1,0/)
res@pmLegendOrthogonalPosF  = -0.20
res@pmLegendParallelPosF    =  0.15

plt = gsn_csm_xy(wks, data(2,:), data(0:1,:), res)

txtres  = True
txtres@txFontHeightF = 0.0145

gsn_text_ndc(wks, Ttext, 0.340, 0.750, txtres)
gsn_text_ndc(wks, Ptext, 0.355, 0.725, txtres)

frame(wks)
