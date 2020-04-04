; ------------------------------------------------ ;
; Created on  Sat Apr 04  18:28:11 2020
; @author : Hyeong-Gyu Kim.
;
; 1. Filename : plot_output.ncl
;
; 2. Variable
;  - Input  : box model result
;  - Output : plot
; ------------------------------------------------ ;

nt = 100
nz = 20

Tout = asciiread("./output.txt", -1, "float")
Tout:= onedtond(Tout, (/nt,nz/))
Tout:= transpose(Tout)
z = fspan(0,10000,nz)

wks = gsn_open_wks("x11", "nonamed")
res = True
res@gsnMaximize = True
; res@vpXF = 0.1
; res@vpYF = 0.8
res@vpWidthF = 0.8
res@vpHeightF = 0.6
res@tmXBMajorOutwardLengthF = 0.0
res@tmXBMinorOutwardLengthF = 0.0
res@tmYLMajorOutwardLengthF = 0.0
res@tmYLMinorOutwardLengthF = 0.0

cres = res
cres@cnFillOn = True
cres@cnLevelSelectionMode = "ExplicitLevels"
cres@cnLevels = ispan(270,295,1)

cres@lbOrientation 	= "vertical"
cres@pmLabelBarHeightF           = 0.6
cres@pmLabelBarWidthF            = 0.06
cres@pmLabelBarOrthogonalPosF    = 0.05
cres@pmLabelBarParallelPosF      = 0.5
cres@lbLabelStride = 2
cres@lbLabelFontHeightF = 0.016
; plt = gsn_csm_contour(wks, Tout, cres)

; make gif plot
wks = gsn_open_wks("pdf", "for_gif")
lres = res
lres@trXMinF = 0.
lres@trXMaxF = 20.
; lres@trYMinF = 
; lres@trYMaxF =

do i = 0, nt-1
    plt = gsn_csm_xy(wks, Tout(:,i), z, lres)
end do    


