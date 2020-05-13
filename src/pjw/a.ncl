radius = asciiread("./log",-1,"float")
radius := reshape(radius, (/476,40,100/))
radius := radius(:,10,:)

r1 = radius(:,0)
r2 = radius(:,10)
r3 = radius(:,20)
r4 = radius(:,30)
r5 = radius(:,40)
r6 = radius(:,50)
r7 = radius(:,60)
r8 = radius(:,70)
r9 = radius(:,80)
r10= radius(:,90)
r11= radius(:,99)
t = ispan(0,475,1)

print(r1)
	;r1 = radius( 0,:)
	;r2 = radius(10,:)
	;r3 = radius(20,:)
	;r4 = radius(30,:)
	;r5 = radius(40,:)
	;r6 = radius(50,:)
	;r7 = radius(60,:)
	;r8 = radius(70,:)
	;r9 = radius(80,:)
	;r10= radius(90,:)
	;r11= radius(99,:)


w = gsn_open_wks("x11","A")
res = True
res@trXMinF = 20
res@trXMaxF = 100
;res@trYMinF
;res@trYMaxF

;res@tiMainString = ""
res@tiXAxisString = "Time [s]"
res@tiYAxisString = "Radius [mm]"

p = gsn_csm_xy(w,t,(/r1,r2,r3,r4,r5,r6,r7,r8,r9,r10,r11/),res)
