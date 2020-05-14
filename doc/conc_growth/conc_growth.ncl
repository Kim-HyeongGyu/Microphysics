radius = asciiread("./log",-1,"float")
printVarSummary(radius)
dt = 42 ;[s]
nt = toint(dimsizes(radius)/100.)
radius := reshape(radius, (/nt,100/))    ; Mass [kg]
; radius := reshape(radius, (/nt,40,100/))
; radius := radius(:,10,:)

bin = ispan(1, 100, 1)

w = gsn_open_wks("pdf","A")
res = True
res@trXMinF = 0
res@trXMaxF = 100
res@trYMinF = min(radius)
res@trYLog  = True

;res@tiMainString = ""
res@tiXAxisString = "Bins"
res@tiYAxisString = "Mass [kg]"

do t = 0, nt-1, 500
    res@gsnRightString = "time: "+(t*dt)+" [s]"
    p = gsn_csm_xy(w,bin,radius(t,:),res)
end do    

