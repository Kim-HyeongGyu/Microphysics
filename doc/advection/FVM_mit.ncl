; This Matlab script solves the one-dimensional convection
; equation using a finite volume algorithm. The
; discretization forward Euler in time.
; URL: http://web.mit.edu/16.90/BackUp/www/pdfs/Chapter16.pdf

; Number of points
Nx = 50
x  = fspan(0,1,Nx+1)
dx = 1./Nx

; Calculate midpoint values of x in each control volume
xmid = 0.5*(x(0:Nx-1) + x(1:Nx))

; Set velocity
u = 1

; Set final time
tfinal = 1.

; Set timestep
CFL = 0.5
dt = CFL*dx/abs(u)

; Set initial condition to U0 = exp(-xˆ2)
; Note: technically, we should average the initial
; distribution in each cell but I chose to just set
; the value of U in each control volume to the midpoint
; value of U0.
U = 0.75*exp(-(((xmid-0.5)/0.1)^2))
t = 0.

wks = gsn_open_wks("x11", "nonamed")
res = True
; Loop until t > tfinal
Ubc = new(dimsizes(U)+2, "float")
do while (t .lt. tfinal)
    Ubc(1:Nx) = U
    Ubc(0) = U(Nx-1)
    Ubc(Nx+1) = U(0)
    ;[U(Nx); U; U(1)]; % This enforces the periodic bc

    ; Calculate the flux at each interface
    F = 0.5* u *( Ubc(1:Nx+1) + Ubc(0:Nx)) - 0.5*abs(u)*( Ubc(1:Nx+1) - Ubc(0:Nx))

    ; Calculate residual in each cell
    R = F(1:Nx) - F(0:Nx-1);

    ; Forward Euler step
    U = U - (dt/dx)*R;

    ; Increment time
    t = t + dt;

    plt = gsn_csm_xy(wks, xmid, U, res)
    print(sum(U)+"  "+sum(Ubc(1:Nx)))
end do

