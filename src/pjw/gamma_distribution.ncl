; Variables
x = ispan(1,100,1)*1.
thetha = 2.0                           ;need to change
k = 1.0                                ;need to change
Nc = 100000000.d
qc = 1.e-6
u = ((10^9)/Nc)
rho = 10^-18                   ;[g/um^3]
pi= 4*atan(1)

; Formula
a = gamma(k)*((thetha)^k) 
b = x^(k-1)
c = exp(-(x/thetha))
GD = (1/a)*b*c                          ;GD = Gamma_Distribution


N0 = (Nc/gamma(u+1))*(rho*(4./3.*pi*Nc)*(u+1)*(u+2)*(u+3)*(1/qc))^((u+1)/3.)
lamda = (rho*(4./3.)*pi*Nc*(u+1)*(u+2)*(u+3)*(1/qc))^(1./3.)

Nr = N0*x^u*exp(-lamda*x)
printMinMax(Nr,False)

; Confirmation of Formula
printMinMax(GD, False)
printMinMax(b, False)
printMinMax(c, False)


; Draw a Graph

w = gsn_open_wks("x11","gamma")

r = True
;r@trYMinF = 0.
;r@trXMinF = 0.
;r@trYMaxF = 0.5
;r@trXMaxF = 20

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
