; Redistribution /*{{{*/ 
undef("redistribution")
function redistribution(m1, m2, n1)
begin
    n2   = n1 * 0.
    nbin = dimsizes(n1)

    do i = 0, nbin-1
        ; Evaporate all drops -> no mass
        if ( m2(i) .lt. m1(0) ) then
            n2(i) = 0
        else
            if ( (i .eq. nbin-1) .and. ( m1(i) .lt. m2(i) ) ) then
                ; Condensation at boundary
                N_i   = m2(i)*n1(i)/m1(i)
                n2(i) = n2(i) + N_i
            else if ( (i .eq. 0) .and. ( m2(i) .lt. m1(i) ) ) then
                ; Evaporation at boundary
                N_i   = m2(i)*n1(i)/m1(i)
                n2(i) = n2(i) + N_i
            else if ( ( m1(i) .lt. m2(i) ) .and. ( m2(i) .lt. m1(i+1) ) ) then
                ; Condensation case
                N_i     = ( m2(i) - m1(i+1) ) \
                        / ( m1(i) - m1(i+1) ) * n1(i)
                N_i1    = n1(i) - N_i
                n2(i)   = n2(i) + N_i
                n2(i+1) = n2(i+1) + N_i1
            else if ( ( m1(i-1) .lt. m2(i) ) .and. ( m2(i) .lt. m1(i) ) ) then
                ; Evaporation case
                N_i     = ( m2(i) - m1(i-1) ) \
                        / ( m1(i) - m1(i-1) ) * n1(i)
                N_i1    = n1(i) - N_i
                n2(i)   = n2(i) + N_i
                n2(i-1) = n2(i-1) + N_i1
            else
                ; No Condensation/Evaporation
                N_i = n1(i)
                N_i1 = 0
                n2(i)   = n2(i) + N_i
            end if    
            end if    
            end if    
            end if    
        end if   
    end do    
    return n2
end ; /*}}}*/

; pprint /*{{{*/
undef("pprint")
procedure pprint(dimension2d[*][*])
; 2-Dimension values print
begin
	option 			= True
	;option@title 	= "row    data"
	option@tspace 	= 2
	;option@row 	 	= True

	prefix = dimsizes(dimension2d(0,:))

	if (typeof(dimension2d) .eq. "integer") then format="I6"   end if
	if (typeof(dimension2d) .eq. "float"  ) then format="F8.3" end if
	if (typeof(dimension2d) .eq. "double" ) then format="F8.3" end if
	
	format = prefix + format

	write_matrix(dimension2d, format, option)
end ; /*}}}*/

; make data
mass1 = (/ 1,   2,         4,    8, 16,  32   /)
num1  = (/ 1,   1,         2,    3,  5,   3   /)
mass2 = (/  1.1,  2.3, 3.9,  7.4,      17, 31 /)    ; grouwh mass (n+1)

num2 = redistribution(mass1, mass2, num1)
data = new((/6,3/), "float")
data(:,0) = mass1
data(:,1) = num1
data(:,2) = num2
print("Mass | Initial num | reassigned num")
pprint(data)
print(sum(num1))
print(sum(num2))
