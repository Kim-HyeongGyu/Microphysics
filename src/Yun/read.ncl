data = fbindirread("./log-normal_dist.gdat",0,40,"double")

wks = gsn_open_wks("x11","res")
res = True
plt = gsn_csm_xy(wks, ispan(1,40,1), data, res)
