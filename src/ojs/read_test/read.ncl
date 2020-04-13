


  f1 = addfile("./data/ERA-Interim_T_2018.nc", "r")
  f2 = addfile("./data/ERA-Interim_Q_2018.nc", "r")


  t = short2flt(f1->t)
  q = short2flt(f2->q)

;  printVarSummary(t)
  print(q(7,:,59,113))
