

data = asciiread("./input_sounding.txt",-1,"float")
nz = 81
data := reshape(data,(/81,3/))
;t = new((/81/),float)
;q = new((/81/),float)

;lev_sfc = data(0,0)
;t_sfc = data(0,0)
lev = data(:,0)
lev(0) = 0
t   = data(:,1)
q   = data(:,2)*0.001


lev!0 = "level"
lev@units = "meter"
lev@long_name = "height"

t!0 = "lev"
t@units = "K"
t@long_name = "Temperature"  
t&lev=lev

q!0 = "lev"
q@units = "kg kg**-1"
q@long_name = "Specific humidity" 
q&lev=lev

;printVarSummary(q)


;exit
;setfileoption("nc","format","netcdf4")
fn1 = "input_sounding_t.nc"
fn2 = "input_sounding_q.nc"
system("rm -f "+fn1); remove if exists
system("rm -f "+fn2); remove if exists
f1 = addfile(fn1,"c")
f2 = addfile(fn2,"c")

;filedimdef(f1,dim_names,dim_sizes,dim_unlimited)
;filedimdef(f2,dim_names,dim_sizes,dim_unlimited)




;f1->lev = lev 
f1->t   = t

;f2->lev = lev 
f2->q   = q
