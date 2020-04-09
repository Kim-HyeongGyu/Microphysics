
;nor = fbindirread("./nor_dis.gdat",0,40,"float")

;print(nor)


pi = 3.141592
nm = 40
ratio = 1.3
n0 = 1000


rb = new(nm+1,"float")
r = new(nm,"float")

rb(0) = 1E-6
std = 2.
ave = 1E-3

do i = 1, nm
    rb(i) = rb(0)*((ratio)^i)
    r(i-1) = (rb(i-1)+rb(i))/2.

end do

print(rb)

nor = new(nm,"float")



do i = 0, nm-1

nor(i)=(n0/((sqrt(2*pi))*log(std)*r(i)))*(exp(-((log(r(i))-log(ave))^2)/2*log(std)^2))

end do


print(nor)

