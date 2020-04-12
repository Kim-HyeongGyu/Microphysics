; Variables

pi    = 4*atan(1)
x     = fspan(1,10000,10000)/1000
N0 = 100
do i = 1,9999
x(i) = x(i-1)*1.00099
end do
print(x(9999))
sigma = 1.5                                 ; need to change
u = 1                                       ; need to change


; Formula
a = x*sigma*sqrt(2.*pi)
b = -((log(x) - u)^2.)
c = 2.*(sigma^2.)
LND = (N0/a) * exp(b/c)                     ; lnd = Log_Normal_Distribution

dx = x

dx(1:) = x(1:) - x(0:9998)

Nc = sum(LND*dx)                            ; Nc = N(r)*dr

rho = 10^-18       ;[g/um^3]

mr = Nc*rho*(4./3.)*pi*(x^3)

qc = sum(mr*dx)

print(qc)


qc2 = rho*(4./3.)*pi*Nc*exp((3./2.)*(2+3*sigma*sigma))
print(qc2)


exit
; Confirmation of formula
printMinMax(LND, False)
printMinMax(b, False)
printMinMax(c, False)

; Draw Graph
w = gsn_open_wks("x11","log")

r = True
r@trYMinF = 0.
r@trYMaxF = 1
r@trXMinF = 0.
r@trXMaxF = 10
r@tmYLMode = "Explicit"
r@tmXTOn = False
r@tmYROn = False
r@gsnFrame = False

r@xyLineColor = "red"
r@tmXMajorGrid = True
r@tmYMajorGrid = True

r@tiMainString = "Log_Normal_Distribution"
r@tiMainFontHeightF = 0.2
r@tiXAxisString = "x"
r@tiYAxisString = "PDF"
r@xyLineThicknesses = 3


p = gsn_csm_xy(w,x,LND,r)


frame(w)
