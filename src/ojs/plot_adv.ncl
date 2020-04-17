

path = "/home/ojs9294/class/mirco_2020/Microphysics/exp/output/"
grid = "stretching"  ; "constant", "stretching"
;ex = "finite_difference"    ;"finite_difference", "finite_volume", "PPM"
ex = "PPM"       ;"finite_difference", "finite_volume", "PPM"
status = "real"  ;"ideal", "real"
vname = "q"        ;"q", "T"
f1 = addfile(path+vname+"_"+status+"_"+grid+"_"+ex+"_out.nc","r")
print(f1)
 dir = status+"_"+grid+"_"+ex
  system ("mkdir "+dir)  
 if (vname .eq. "q") then
 var = f1->q
 else
 var = fi->T
 end if
 z = f1->hgt
 t = f1->time
 nt = dimsizes(t)
; print(nt)
; printVarSummary(q)


  out = "png"
 do tt = 0,nt-1

   outfile = dir+"/"+vname+"_"+status+"_"+ex+tt
   wks = gsn_open_wks(out,outfile)
   res = True
  res@gsnMaximize = True
  res@vpHeightF= 0.4                    ; change aspect ratio of plot
  res@vpWidthF = 0.8   

  res@tiMainString = ex+" advection"
  res@gsnLeftString     = "t = "+ tt
  
    if (vname .eq. "q") then
        res@tiYAxisString = "q(kg/kg)"
       if (status .eq. "ideal") then
           res@trYMinF = -0.1
           res@trYMaxF = 0.5
       else
           res@trYMinF = -0.02
           res@trYMaxF = 0.14
       end if
    else
        res@tiYAxisString = "T(K)"
           res@trYMinF = 200
           res@trYMaxF = 300
    end if

  res@tiXAxisString = "z"


  plot = gsn_csm_xy(wks, z, var(tt,:),res)

 end do



