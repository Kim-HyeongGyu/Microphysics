; Log normal distribution PDF /*{{{*/
undef("log_normal_distribution")
function log_normal_distribution(x, Nc, qc)
begin
    rho = 1.e3                  ;[kg m-3]
    pi  = 4*atan(1)

    N0  = Nc
    r0  = ( (3./(4.*pi))*(qc/Nc)*(1./rho) )^(1./3.)
    u   = log(r0)
    sigma = 2.5
    ; sigma = sqrt(((2./3.)*log((3./4.)/pi/rho/N0)-2*u)/3.)
    Nr  = N0 / (sqrt(2*pi)*sigma*x) \
        * exp(-( (log(x)-u)^2 / (2*sigma^2) ))

    CDF = 0.5 + 0.5*erf((log(x)-u)/(sqrt(2)*sigma))
    Nr@CDF = CDF(dimsizes(CDF)-1)
    return Nr
end ; /*}}}*/

; read data
data = asciiread("./Nr.txt", 80, "double")
x   = data(:39)
; x = ispan(1,1000000,1)*1.e-6
Nr2 = data(40:)

; Variables
Nc  = 1.e3
qc  = 1.

; Formula
Nr = log_normal_distribution(x, Nc, qc)
; print("NCL      Fortran     Difference")
; print(Nr+"    "+Nr2+"    "+(Nr-Nr2))
; exit

; verify PDF
dx  = x(1:) - x(:dimsizes(x)-2)
dx := array_append_record(x(0), dx, 0)
cdf = Nr*dx
print(sum(cdf))
print(Nr@CDF)


w = gsn_open_wks("x11","gamma")

r = True
r@trXMinF = 1e-6
r@trXMaxF = 10
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
