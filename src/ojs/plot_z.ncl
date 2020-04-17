

path = "/home/ojs9294/class/mirco_2020/Microphysics/exp/output/"
ex = "ppm"
f1 = addfile(path + "T_"+ex+"_out.nc","r")
f2 = addfile(path + "T_constant_"+ex+"_out.nc","r")

  

 q = f1->T
 z1 = f1->hgt
 z2 = f2->hgt
 t = f1->time
 nz = dimsizes(z1)
 data = new((/2,nz/),float)

 data(0,:) = z1
 data(1,:) = z2

  out = "png"

;   outfile = ex+"/"+ex+tt
   wks = gsn_open_wks(out,"zz")
   res = True
  res@gsnMaximize = True
  res@vpHeightF= 0.4                    ; change aspect ratio of plot
  res@vpWidthF = 0.8   

;  res@tiMainString = ex+" advection"
;  res@gsnLeftString     = "t = "+ tt
   xz = ispan(1,20,1)

;  res@trYMinF = 200
;  res@trYMaxF = 300
;  res@trYMinF = -3.0
;  res@trYMaxF = 3.0

;  res@tiYAxisString = "q(kg/kg)"
; res@tiYAxisString = "T(K)"
; res@tiXAxisString = "z"


  plot = gsn_csm_xy(wks, xz, data,res)




