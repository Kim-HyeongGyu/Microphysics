; ------------------------------------------------ ;
; Created on  Mon Jun 01  15:13:38 2020
; @author : Hyeong-Gyu Kim.
;
; 1. Filename : bin2nc.ncl
;
; 2. Variable
;  - Input  : Model output (Binary format)
;  - Output : Netcdf format file
;
; 3. Description
;  - Convert from binary to NetCDF
; ------------------------------------------------ ;

; from nml_info.nml
nt   = 5
nz   = 40
nbin = 100
fout = "processed_model.nc"

;===================== DO not touch below code ====================
print("Converting binary file to NetCDF ...")
time            = fbinrecread(            "./time.bin", 0,     nt, "double")
z_full          = fbinrecread(          "./z_full.bin", 0,     nz, "double")
z_half          = fbinrecread(          "./z_half.bin", 0,   nz+1, "double")
radius          = fbinrecread(          "./radius.bin", 0,   nbin, "double")
radius_boundary = fbinrecread( "./radius_boundary.bin", 0, nbin+1, "double")
mass            = fbinrecread(            "./mass.bin", 0, (/nz,  nbin/), "double")
mass_boundary   = fbinrecread(   "./mass_boundary.bin", 0, (/nz,nbin+1/), "double")

Prs    = new( (/nt,nz          /), "double" )     ; [time] x [lev]
T      = new( (/nt,nz          /), "double" )     ; [time] x [lev]
THETA  = new( (/nt,nz          /), "double" )     ; [time] x [lev]
qv     = new( (/nt,nz          /), "double" )     ; [time] x [lev]
W      = new( (/nt,nz+1        /), "double" )     ; [time] x [lev+1]
Nr     = new( (/nt,nz  , nbin  /), "double" )     ; [time] x [lev] x [nbin]
dm_dt  = new( (/nt,nz  , nbin  /), "double" )     ; [time] x [lev] x [nbin]
dmb_dt = new( (/nt,nz  , nbin+1/), "double" )     ; [time] x [lev] x [nbin+1]
do nrec = 0, nt-1
         Prs(nrec,:) = fbinrecread(   "./Prs.bin", nrec, (/nz       /), "double")
           T(nrec,:) = fbinrecread(     "./T.bin", nrec, (/nz       /), "double")
       THETA(nrec,:) = fbinrecread( "./THETA.bin", nrec, (/nz       /), "double")
          qv(nrec,:) = fbinrecread(    "./qv.bin", nrec, (/nz       /), "double")
           W(nrec,:) = fbinrecread(     "./W.bin", nrec, (/nz+1     /), "double")
        Nr(nrec,:,:) = fbinrecread(    "./Nr.bin", nrec, (/nz,nbin  /), "double")
     dm_dt(nrec,:,:) = fbinrecread( "./dm_dt.bin", nrec, (/nz,nbin  /), "double")
    dmb_dt(nrec,:,:) = fbinrecread("./dmb_dt.bin", nrec, (/nz,nbin+1/), "double")
end do    

; Build coordinate variables
time!0 = "time"
time@long_name = "time (index x dt)"
time@units = "s"

z_full!0 = "z_full"
z_full@long_name = "Full coordinate in model (nz)"
z_full@units = "m"

z_half!0 = "z_half"
z_half@long_name = "Half coordinate in model (nz+1)"
z_half@units = "m"

bin = ispan(1, nbin, 1)  
bin_b = ispan(0, nbin, 1)

radius!0 = "radius"
radius@long_name = "droplet radius"
radius@units = "m"

radius_boundary!0 = "radius_boundary"
radius_boundary@long_name = "droplet radius at boundary"
radius_boundary@units = "m"

mass!0 = "radius"
mass!1 = "bin"
mass@long_name = "droplet radius"
mass@units = "m"

mass_boundary!0 = "radius_boundary"
mass_boundary!1 = "bin_b"
mass_boundary@long_name = "droplet radius at boundary"
mass_boundary@units = "m"

Prs!0 = "time"
Prs!1 = "z_full"
Prs&time = time
Prs&z_full = z_full
Prs@long_name = "Pressure"
Prs@units = "hPa"

copy_VarCoords(Prs, T)
T@long_name = "Temperature"
T@units = "K"

copy_VarCoords(Prs, THETA)
THETA@long_name = "Potential temperature"
THETA@units = "K"

copy_VarCoords(Prs, qv)
qv@long_name = "Mixing ratio"
qv@units = "kg kg-1"

copy_VarCoords(Prs(:,0), W(:,0))
W!1 = "z_half"
W&z_half = z_half
W@long_name = "Vertical wind speed"
W@units = "m s-1"

Nr!0 = "time"
Nr!1 = "z_full"
Nr!2 = "bin"
Nr&time = time
Nr&z_full = z_full
Nr&bin    = bin
Nr@long_name = "Number of droplets"
Nr@units = "# m-3 (Maybe?)"

copy_VarCoords(Nr, dm_dt)
dm_dt@long_name = "Tendency of mass"
dm_dt@units = "kg s-1"

copy_VarCoords(Nr(:,:,0), dmb_dt(:,:,0))
dmb_dt!2 = "bin_b"
dmb_dt&bin_b = bin_b
dmb_dt@long_name = "Tendency of mass"
dmb_dt@units = "kg s-1"

system("rm -f "+fout)
setfileoption("nc","Format","LargeFile")    ; > 2GB
fc = addfile(fout, "c")
; fc->time            = time
; fc->z_full          = z_full
; fc->z_half          = z_half
; fc->radius          = radius
; fc->radius_boundary = radius_boundary
fc->mass            = mass
fc->mass_boundary   = mass_boundary   
fc->Prs             = Prs             
fc->T               = T               
fc->THETA           = THETA           
fc->qv              = qv              
fc->W               = W               
fc->Nr              = Nr              
fc->dm_dt           = dm_dt           
fc->dmb_dt          = dmb_dt          

print("Done! Output --> "+fout)
