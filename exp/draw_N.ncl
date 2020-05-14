data1 = asciiread("./rb.txt", -1, "float")   ; rb
data2 = asciiread("./Nr.txt", -1, "float")   ; Nr
;printVarSummary(data)
nbin = 100
ntime = 23808
rb = reshape(data1, (/nbin+1/))
Nr = reshape(data2, (/ntime,nbin/))

dlnr= new((/nbin/),float)
r= new((/nbin/),float)
do i = 0, nbin-1
    r(i) = (rb(i+1)+rb(i))/2
    dlnr(i) = log(rb(i+1)/rb(i))
end do



;===========================plot

do i = 0, ntime,100
wks = gsn_open_wks("png", "./Nr_100/Nr_"+i)
res = True
;res@trXMinF  = rb(1)
;res@trXMaxF  = rb(nbin)
res@trYMinF  = 0.
res@trYMaxF  = 2.e+8
res@trXLog = True

res@tmXTOn = False
res@tmYROn = False
res@xyLineThicknessF = 3.0
;res@xyXStyle = "Log"
res@tiXAxisString = "r[m]"
res@tiYAxisString = "dN/dlnr[#/m]"

res@gsnLeftString = "t = "+ i





;res@tmXBMode = "Manual"
;res@tmXBMode = "Explicit"

;res@tmXBTickStartF= rb(1)
;res@tmXBTickEndF = rb(nbin)
;res@tmXBTickSpacingF = 1.e-03
;res@tmXBLabels = rr 
;do i = 0, ntime, 10
    plt = gsn_csm_xy(wks, r, Nr(i,:)/dlnr(:), res)
;    plt = gsn_csm_xy(wks, r(:), Nr(i,:)/dlnr(:), res)
end do    

