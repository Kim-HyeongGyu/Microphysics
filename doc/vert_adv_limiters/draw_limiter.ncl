data1 = asciiread("./Nr_from_Colella.txt", -1,"float")
data2 = asciiread("./Nr_from_Lin.txt", -1,"float")
data1 := reshape(data1,(/119999,40/))
data2 := reshape(data2,(/119999,40/))

wks = gsn_open_wks("png","./Nr_limiters/Nr_limiters")
res = True

res@gsnFrame    = False
res@gsnMaximize = True
res@vpHeightF   = 0.5
res@vpWidthF    = 0.5

res@xyLineColors        = (/"black","red"/)
res@xyDashPatterns      = (/0.,2./)
res@xyLineThicknesses   = (/4.,4./)

;res@tiMainString    = "Limiters from Colella Woodward (1984)"
res@tiXAxisString   = "z [m]"
res@tiYAxisString   = "N(r)"

res@trYMinF = 0.
res@trYMaxF = 7.e6
res@trXMinF = 0.
res@trXMaxF = 20. 

res@tmXTOn  = False
res@tmYROn  = False
res@tmXBMode    = "Explicit"
res@tmXBValues  = (/0, 10, 20/)
res@tmXBLabels  = (/"0","500","1000"/)
res@tmYLMode    = "Explicit"
res@tmYLValues  = (/0, 1e6, 2e6, 3e6, 4e6, 5e6, 6e6, 7e6/)
res@tmYLLabels  = (/"0.0~F34~4~F21~10~S~6","1.0~F34~4~F21~10~S~6","2.0~F34~4~F21~10~S~6","3.0~F34~4~F21~10~S~6","4.0~F34~4~F21~10~S~6","5.0~F34~4~F21~10~S~6","6.0~F34~4~F21~10~S~6","7.0~F34~4~F21~10~S~6"/) 

res@pmLegendDisplayMode = "Always"
res@pmLegendSide        = "Top"
res@pmLegendHeightF     = 0.04
res@pmLegendWidthF      = 0.15
res@lgPerimOn           = False
res@lgLabelFontHeightF  = 0.01
res@xyExplicitLegendLabels  = (/"  Colella and Woodward (1984)","  Lin (2003)"/)
res@lgItemOrder         = (/1,0/)
res@pmLegendOrthogonalPosF  = -0.11
res@pmLegendParallelPosF    = 0.3

do i = 0, 119999, 1000
plt = gsn_csm_xy(wks, ispan(1,40,1), (/data1(i,:), data2(i,:)/), res)

txtres  = True
txtres@txFontHeightF = 0.02
gsn_text_ndc(wks, "t = "+i+" s", 0.85, 0.88, txtres)
frame(wks)
end do
