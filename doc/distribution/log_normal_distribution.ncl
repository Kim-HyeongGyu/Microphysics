; Gamma distribution PDF /*{{{*/
undef("log_normal_distribution")
function log_normal_distribution(rb, Nc, qc)
begin
    rho = 1.e3                  ;[kg m-3]
    pi  = 4*atan(1)

    r   = ( rb(1:)+rb(:dimsizes(rb)-2) ) / 2.d
    D   = 2*r
    dr  = rb(1:) - rb(:dimsizes(rb)-2)
    dD  = 2*dr

    N0  = Nc
    r0_max = ( qc / (Nc*rho*(4./3.d)*pi) )^(1./3.)
    r0     = 1.e-5  ; r0 < r0_max
    ; r0     = r0_max  ; r0 < r0_max
    if (r0 .gt. r0_max) then
        print("Note! r0 > r0_max")
        exit
    end if    

    u     = log(r0)
    sigma = sqrt((2./9.)*log(qc/(rho*(4./3.)*pi*(r0^3.)*Nc)))

    Nr  = N0/(sqrt(2*pi)*r*sigma)   \
        * exp(- ( (log(r)-u)^2 )    \
              / (    2*sigma^2 ) ) * dr
    m   = (4./3.)*pi*rho*(r^3)*Nr

    dN_dlnr = Nr
    do i = 0, dimsizes(r)-1
        dN_dlnr(i) = Nr(i)/log(rb(i+1)/rb(i))
    end do    

    ; CDF  = Nrdr
    ; do i = 0, dimsizes(x)-1
    ;     CDF(i) = (1./gamma(u+1))*gammainc(lamda*x(i),u+1) 
    ; end do    
    ; Nrdr@CDF = CDF(dimsizes(CDF)-1)
    dN_dlnr@r  = r
    dN_dlnr@Nc = sum(Nr)
    dN_dlnr@qc = sum(m)

    return dN_dlnr
end ; /*}}}*/

; Variables
Nc   = 1.e8  ; [m-3]
rmin = 1.e-6 ; [m]
rmax = 1.e-2 ; [m]
nbin = 120   ; [#]
rratio = (rmax/rmin)^(1./nbin)  ; = 1.096478

; radius
rb = new(nbin+1, "float")
do i = 0, nbin
    rb(i) = rmin*(rratio^i)
end do    

; Formula
qc  = 0.5e-3   ; [kg kg-1]
Nr1 = log_normal_distribution(rb, Nc, qc)
qc  = 1.e-3
Nr2 = log_normal_distribution(rb, Nc, qc)
qc  = 2.e-3
Nr3 = log_normal_distribution(rb, Nc, qc)

; verify PDF
print("   Nc   |   qc")
print(Nr1@Nc+"  "+Nr1@qc)
print(Nr2@Nc+"  "+Nr2@qc)
print(Nr3@Nc+"  "+Nr3@qc)

wks = gsn_open_wks("x11","gamma")

res = True
res@vpWidthF = 0.8
res@vpHeightF = 0.4
res@gsnMaximize = True

res@trXMinF = 1.e-6
res@trXMaxF = 0.02
res@trYMinF = 0.
; res@trYMaxF = 2.
res@trXLog  = True

; res@tmXBMode   = "Explicit"
; res@tmXBValues = ispan(0, 20, 5)*1.e-3
; res@tmXBLabels = res@tmXBValues
; res@tmYLMode   = "Explicit"
; res@tmYLValues = ispan(0, 35, 5)*1.e7
; res@tmYLLabels = res@tmYLValues
res@tmXTOn = False
res@tmYROn = False

res@xyLineColors = (/"red", "green", "blue"/)
res@xyLineThicknesses = (/5, 5, 5/)
res@xyDashPatterns = (/0, 0, 0/)
res@tmXMajorGrid                = True      ; X-axis grid on
res@tmYMajorGrid                = True      ; Y-axis grid on
res@tmXMajorGridLineDashPattern = 2         ; X-axis grid style
res@tmYMajorGridLineDashPattern = 2         ; Y-axis grid style
res@tmXMajorGridThicknessF      = 0.8       ; X-axis grid thick
res@tmYMajorGridThicknessF      = 0.8       ; Y-axis grid thick      
res@tmXMajorGridLineColor		= "grey40"  ; X-axis grid color
res@tmYMajorGridLineColor		= "grey40"  ; Y-axis grid color

res@gsnLeftString = "Log-normal Distribution"
res@gsnLeftStringFontHeightF = 0.02
res@tiXAxisString = "r [m]"
res@tiYAxisString = "N(r) x dr"
plt = gsn_csm_xy(wks,Nr1@r,(/Nr1, Nr2, Nr3/),res)
; p = gsn_csm_xy(w,x,Nr@CDF,r)
