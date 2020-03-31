; Tout = fbinread("./output.dat", (/100,20/), "float")
; Tout = fbinrecread("./output.dat", 0, 20, "float")
; Tout = fbindirread("./output.dat", 0, (/100,20/), "float")
Tout = asciiread("./output.txt", -1, "float")
Tout:= reshape(Tout, (/100,20/))
printVarSummary(Tout)
; print(Tout(:2,:))
; printMinMax(Tout(:2,:), False)
; exit

Tout := transpose(Tout)
ztop = 16.*1000.
zbot = 0
dz = (ztop-zbot) / 20.

z = new(20, "float")
z(0) = zbot
do i = 1, 19
    z(i) = z(i-1) + dz
end do    

wks = gsn_open_wks("x11", "nonamed")
res = True
res@cnFillOn = True
res@cnLevelSelectionMode = "ExplicitLevels"
res@cnLevels = ispan(270,295,1)
plt = gsn_csm_contour(wks, Tout, res)
