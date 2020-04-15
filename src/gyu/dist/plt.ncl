; Variables
x   = ispan(1,100000,1)*1.e-3
Nc  = 1.e8
qc  = 1.e-6
rho = 10^-18                   ;[g/um^3]
pi  = 4*atan(1)

; Formula
u     = ((10^9)/Nc)
umul  = (u+1)*(u+2)*(u+3)
N0    = (Nc/gamma(u+1)) \ 
      * ( (rho/qc)*(4./3.)*pi*Nc*umul)^((u+1)/3.)
lamda = ( (rho/qc)*(4./3.)*pi*Nc*umul)^(1./3.)

Nr = N0*x^u*exp(-lamda*x)
Nr2 = asciiread("./Nr.txt", 100000, "double")
; print("NCL      Fortran     Difference")
; print(Nr+"    "+Nr2+"    "+(Nr-Nr2))
; exit

w = gsn_open_wks("x11","gamma")

r = True
; r@trXMinF = 0.001
; r@trXMaxF = 10
; r@trYMinF = 0.
; ; r@trYMaxF = 2.
; r@trXLog  = True

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

p = gsn_csm_xy(w,x,Nr2,r)

frame(w)
