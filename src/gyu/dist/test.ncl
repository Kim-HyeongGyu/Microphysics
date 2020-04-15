; Variables
data = asciiread("./Nr.txt", 80, "double")
x = data(:39)
Nr = data(40:)
; print("NCL      Fortran     Difference")
; print(Nr+"    "+Nr2+"    "+(Nr-Nr2))
; exit
printMinMax(Nr, False)

w = gsn_open_wks("x11","gamma")

r = True
r@trXMinF = 1.e-6
r@trXMaxF = 1.
r@trYMinF = 0.
; r@trYMaxF = 2.
r@trXLog  = True

r@tmYLMode = "Explicit"
r@tmXTOn = False
r@tmYROn = False
r@gsnFrame = False

r@xyLineColor = "red"
r@tmXMajorGrid = True
r@tmYMajorGrid = True

r@tiMainString = "Gamma Distribution"
r@tiMainFontHeightF = 0.2
r@xyLineThicknesses = 3

p = gsn_csm_xy(w,x,Nr,r)

frame(w)
