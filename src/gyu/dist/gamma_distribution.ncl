; Gamma distribution PDF /*{{{*/
undef("gamma_distribution")
function gamma_distribution(x, Nc, qc)
begin
    rho = 1.e3                  ;[kg m-3]
    pi  = 4*atan(1)

    u     = 1.
    ; u     = ((10^9)/Nc)
    umul  = (u+1)*(u+2)*(u+3)
    N0    = (Nc/gamma(u+1)) \ 
          * ( (rho/qc)*(4./3.)*pi*Nc*umul)^((u+1)/3.)
    lamda = ( (rho/qc)*(4./3.)*pi*Nc*umul)^(1./3.)
    Nr    = N0*x^u*exp(-lamda*x)

    CDF  = Nr
    do i = 0, dimsizes(x)-1
        CDF(i) = (1./gamma(u+1))*gammainc(lamda*x(i),u+1) 
    end do    
    Nr@CDF = CDF

    return Nr
end ; /*}}}*/

; read data
data = asciiread("./Nr.txt", 80, "double")
x   = data(:39)
; x   = ispan(1, 100000, 1)*1.e-4
Nr2 = data(40:)

; Variables
Nc  = 1.e3
qc  = 1e-6

; Formula
Nr = gamma_distribution(x, Nc, qc)
; print("NCL      Fortran     Difference")
; print(Nr+"    "+Nr2+"    "+(Nr-Nr2))
; exit

; verify PDF
dx  = x(1:) - x(:dimsizes(x)-2)
dx := array_append_record(x(0), dx, 0)
cdf = Nr*dx
; print(sum(cdf))
print(max(Nr@CDF))

w = gsn_open_wks("x11","gamma")

r = True
r@trXMinF = 1.e-6
r@trXMaxF = 10
r@trYMinF = 0.
; r@trYMaxF = 2.
r@trXLog  = True

r@tmYLMode = "Explicit"
r@tmXTOn = False
r@tmYROn = False

r@xyLineColor = "red"
r@tmXMajorGrid = True
r@tmYMajorGrid = True

r@tiMainString = "Gamma Distribution"
r@tiMainFontHeightF = 0.2
r@xyLineThicknesses = 3

p = gsn_csm_xy(w,x,Nr,r)
p = gsn_csm_xy(w,x,Nr@CDF,r)
