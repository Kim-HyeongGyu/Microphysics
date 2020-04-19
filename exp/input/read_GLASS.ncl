f = addfile("GLASS.nc", "r")
z = f->alt  ; Altitude     [m]
p = f->pres ; Pressure     [hPa]
T = f->tdry ; Temperature  ['C]
w = f->mr   ; Mixing ratio [g/kg]
