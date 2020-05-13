undef("terminal_velocity") ; /*{{{*/
function terminal_velocity(radius, T, P)
; Input
; - radius : [um]
; - T      : Temperature [K]
; - P      : Pressure [hPa]
;
; Reference 
; - Beard (1976)
begin
    
    d0 = radius * 2   ; radius -> diameter [um]
    _FillValue = 0.

    ; constant
    R   = 287.       ; [J kg-1 K-1]
    g   = 9.8        ; [m s-2]
    T0  = 293.15     ; [K]
    P0  = 1013.25    ; [hPa]
    l0  = 6.62e-6    ; [cm]
    mu0 = 0.0001818  ; [g cm-1 s-1]

    rho_liquid = 1000.      ; [kg m-3] water density
    rho_air    = (P*100.)/(R*T)    ; [kg m-3] air density  ; p = rho R T -> rho = P/RT
    drho = rho_liquid - rho_air    ; drop - air

    ; dynamic viscosity ( Approximate formula, See Yau (1996) - 102p )
    mu = 1.72e-5 * ( 393/(T+120.) ) * ( (T/273)^(3./2.) )   ; [kg m-1 s-1]

    l   = l0 * (mu/mu0) * (P0/P) * sqrt(T/T0)  ; [cm] mean free path of air molecules
    C1  = drho * g / (18*mu)        ; [m-1 s-1] = [kg m-3] * [m s-2] / [kg m-1 s-1]
    Csc = 1 + 2.51*l/d0             ; [dimensionless]
    ; Calculate terminal velocity in each regime
    if (d0 .lt. 0.5 ) then
        Vt = _FillValue     ; ignore
    else if (d0 .le. 19) then
        ; Regime 1
        Vt = C1 * Csc * (d0*1.e-6)^2   ; [m s-1] = [m-1 s-1] [dimensionless] [um^2]
    else if (d0 .le. 1.07e3) then
        ; Regime 2
        b0  = -0.318657e+1
        b1  =  0.992696
        b2  = -0.153193e-2
        b3  = -0.987059e-3
        b4  = -0.578878e-3
        b5  =  0.855176e-4
        b6  = -0.327815e-5
        C2  = 4 * rho_air * drho* g / ( 3 * mu^2 )  ; [] = [kg2 m-6] [m s-2] / [um2]
        Da  = C2 * (d0*1e-6)^3  ; Davies number [kg2 s-2] = [kg2 m-3 s-2] [m-3]
        X   = log( Da )
        Y   = b0 + b1*X + b2*X^2 + b3*X^3 + b4*X^4 + b5*X^5 + b6*X^6
        Re  = Csc * exp(Y)              ; Reynolds number
        Vt  = mu * Re / (rho_air * d0*1.e-6)
        ; print("TODO: Regime 2 problem, Vt == 0")
        ; exit
    else if (d0 .le. 7.e3) then
        ; Regime 3
        b0  = -0.500015e+1
        b1  =  0.523778e+1
        b2  = -0.204914e+1
        b3  =  0.475294
        b4  = -0.542819e-1
        b5  =  0.238449e-2

        ; surface tension [N m-1]
        ; sigma = 7.5 * 1e-2  ; Yau (1996) - 85p
        ; See Yau (1996) - 6.9 problem
        Cl = -1.55 * 1e-4   ; [N m-1 K-1]
        C2 = 0.118          ; [N m-1]
        sigma = Cl*T + C2   ; Resonable at -20 ~ 20 [K] temperature

        C3  = 4 * drho * g / (3*sigma)
        Bo  = C3 * (d0*1.e-6)^2.       ; modified Bond number
        Np  = sigma^3. * rho_air^2. / (mu^4. * drho * g)
        X   = log( Bo * Np^(1./6.) )
        Y   = b0 + b1*X + b2*X^2 + b3*X^3 + b4*X^4 + b5*X^5

        Re  = Np^(1./6.) * exp(Y)
        Vt  = mu * Re / (rho_air*(d0*1.e-6))
    else
	d0  = 7.e3
        b0  = -0.500015e+1
        b1  =  0.523778e+1
        b2  = -0.204914e+1
        b3  =  0.475294
        b4  = -0.542819e-1
        b5  =  0.238449e-2

        ; surface tension [N m-1]
        ; sigma = 7.5 * 1e-2  ; Yau (1996) - 85p
        ; See Yau (1996) - 6.9 problem
        Cl = -1.55 * 1e-4   ; [N m-1 K-1]
        C2 = 0.118          ; [N m-1]
        sigma = Cl*T + C2   ; Resonable at -20 ~ 20 [K] temperature

        C3  = 4 * drho * g / (3*sigma)
        Bo  = C3 * (d0*1.e-6)^2.       ; modified Bond number
        Np  = sigma^3. * rho_air^2. / (mu^4. * drho * g)
        X   = log( Bo * Np^(1./6.) )
        Y   = b0 + b1*X + b2*X^2 + b3*X^3 + b4*X^4 + b5*X^5

        Re  = Np^(1./6.) * exp(Y)
        Vt  = mu * Re / (rho_air*(d0*1.e-6))
    end if    
    end if    
    end if    
    end if    
     
    return Vt
end ; /*}}}*/

radius = ispan(10, 5000, 10)*1. ; [um]
; radius = 5.e2     ; [um]
temp1   = 293.15    ; [K]
prs1    = 1013    ; [hPa]

temp2   = 293.15
prs2    = 900

temp3   = 288.15
prs3    = 800

temp4   = 283.15
prs4    = 700

temp5   = 273.15
prs5    = 600

temp6   = 263.15
prs6    = 500

Vt1   = radius
nbin = dimsizes(radius)
do r = 0, nbin-1
    Vt1(r) = terminal_velocity(radius(r), temp1, prs1)
end do   

Vt2  = radius
nbin = dimsizes(radius)
do r = 0, nbin-1
   Vt2(r) = terminal_velocity(radius(r), temp2, prs2)
end do

Vt3  = radius
nbin = dimsizes(radius)
do r = 0, nbin-1
   Vt3(r) = terminal_velocity(radius(r), temp3, prs3)
end do

Vt4  = radius
nbin = dimsizes(radius)
do r = 0, nbin-1
   Vt4(r) = terminal_velocity(radius(r), temp4, prs4)
end do


Vt5  = radius
nbin = dimsizes(radius)
do r = 0, nbin-1
   Vt5(r) = terminal_velocity(radius(r), temp5, prs5)
end do


Vt6  = radius
nbin = dimsizes(radius)
do r = 0, nbin-1
   Vt6(r) = terminal_velocity(radius(r), temp6, prs6)
end do


wks = gsn_open_wks("x11", "nonamed")
res = True

res@tiMainString  = "Terminal velocity of raindrops at pressure"
res@tiXAxisString = "Diameter [mm]"
res@tiYAxisString = "Terminal Velocity [m s-1]"

res@trXMinF = 0.
res@trXMaxF = 9.
res@trYMinF = 0.
res@trYMaxF = 13.

res@xyLineThicknessF = 3


diameter = 2 * radius * 1.e-3   ; [um] -> [mm]


p = gsn_csm_xy (wks, diameter, (/Vt1, Vt2, Vt3, Vt4, Vt5, Vt6/), res)

