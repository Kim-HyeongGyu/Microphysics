      implicit double precision (a-h,o-z)
      parameter (nkt=120,nka=120,nb=20,n=150,nktt=80,nkaa=40)
      common /cb29/ c(nkt,nkt),elnr(nkt),ima(nkt,nkt),jmin
      common /cb50/ enw(nka),ew(nkt),rn(nka),rw(nka,nkt),en(nka),
     &              e(nkt),dew(nkt),rq(nka,nkt),kg,kgp,kd,kdp
      common /cb52/ fpi(nka),fpj(nka),f(nka,nkt),g(nkt),fs(nkt)
      common /cb79/ xm1pl,dlnr,rl(nkt),ax,axn,i4
      common/courn/cn(nka,nka),icn(nka,nka)
      common/kern/ck(nkt,nkt)
      common /co2d/aerm,wini,fogtype,filename
      real if1,if2,ixm2,irh,izk,ix0,ix1,ix2
      dimension if1(nkaa,nktt),if2(nkaa,nktt)
      character *1 ctime
      data dt,ilev,lday1,lst1,lmin1 /10.,3,0,6,0/
c coalbo.dat produced by stratus model as file ci .out
      open (04,file='coalbo.dat',status='old',form='unformatted')
c f2din.dat: input 2d particle spectrum
      open (05,file='f2din.dat',status='old',form='unformatted')
c f2dbo.out: output 2d particle spectrum
      open (06,file='f2dbo.out',status='new',form='unformatted')
c coldat.out: number of collisions per timestep
      open (07,file='coldat.out',status='new')
      read (04) cn,ck,c,elnr,rn,en,e,rl,dlnr,ax,axn,ima,jmin,icn
      close (04)
      print *,'enter timestep'
      read *,dt
      print *,'enter day, hour, minute'
      read *,lday1,lst1,lmin1
      print *,'enter level'
      read *,ilev
      print *,'timestep: ',dt
      print *,'times: ',lday1,lst1,lmin1
      print *,'level: ', ilev
      do idt=1,1000
         if (ilev.eq.1) then
            read (5,end=2000) lday,lst,lmin,kf,ixm2,irh,izk,if2
            read (5,end=2000) lday,lst,lmin,kf,ix0,ix1,ix2,if1
            read (5,end=2000) lday,lst,lmin,kf,ix0,ix1,ix2,if1
         elseif (ilev.eq.2) then
            read (5,end=2000) lday,lst,lmin,kf,ix0,ix1,ix2,if1
            read (5,end=2000) lday,lst,lmin,kf,ixm2,irh,izk,if2
            read (5,end=2000) lday,lst,lmin,kf,ix0,ix1,ix2,if1
         elseif (ilev.eq.3) then
            read (5,end=2000) lday,lst,lmin,kf,ix0,ix1,ix2,if1
            read (5,end=2000) lday,lst,lmin,kf,ix0,ix1,ix2,if1
            read (5,end=2000) lday,lst,lmin,kf,ixm2,irh,izk,if2
         endif
         if (lday.eq.lday1.and.lst.eq.lst1.and.lmin.eq.lmin1) then
            print *,lday,lst,lmin,kf,ixm2,irh,izk
            aerm=0.
            wini=0.
            do jt=1,nkt
            do ia=1,nka
               f(ia,jt)=0.
            enddo
            enddo
            do jt=1,nktt
               do ia=1,nkaa
                  f(ia,jt)=if2(ia,jt)
               enddo
               do ia=1,nkaa
                  aerm=aerm+f(ia,jt)*en(ia)*1.d09
                  wini=wini+f(ia,jt)*e(jt)
               enddo
            enddo
            print *,'initial water and aerosol mass',wini*1000.,aerm
            go to 2000
         endif
      enddo
 2000 continue
      close (05)
      call cofm (dt)
      close (06)
      close (07)
      stop
      end

      subroutine cofm (dt)
      implicit double precision (a-h,o-z)
      parameter (nkt=120,nka=120,nb=20,n=150)
      common /cb29/ c(nkt,nkt),elnr(nkt),ima(nkt,nkt),jmin
      common /cb50/ enw(nka),ew(nkt),rn(nka),rw(nka,nkt),en(nka),
     &              e(nkt),dew(nkt),rq(nka,nkt),kg,kgp,kd,kdp
      common /cb52/ fpi(nka),fpj(nka),f(nka,nkt),g(nkt),fs(nkt)
      common /cb79/ xm1pl,dlnr,rl(nkt),ax,axn,i4
      common/kern/ck(nkt,nkt)
      common/courn/cn(nka,nka),icn(nka,nka)
      common /datco/ coldat(70,10),icol
      common /co2d/aerm,wini,fogtype,filename
      real if2,ixm2,irh,izk
      dimension if2(nka,nkt)
      character*7 filename
      character *1 fogtype
      data tmax/3600./
c      call fallg
c      call trkern(0)
c begin one-dimensional coalescence
c      do jt=1,nkt
c      do ia=2,nka
c         f(1,jt)=f(1,jt)+f(ia,jt)
c         f(ia,jt)=0.
c      enddo
c      enddo
c end one-dimensional coalescence
      do jt=1,nkt
         i4=jt
         if (rl(jt).gt.1.) go to 2000
      enddo
 2000 continue
      lday=0
      lst=0
      lmin=0
      do ia=1,nka
      do jt=1,nkt
         if2(ia,jt)=f(ia,jt)
      enddo
      enddo
      write (6) lday,lst,lmin,kf,ixm2,irh,izk,if2
 6100 format (5e16.8)
      do ij=1,70
      do ilk=1,10
         coldat(ij,ilk)=0.0
      enddo
      enddo
      nt=int(tmax/dt)
      tlmin=1.d-5
      do ij=1,nt
      icol=lmin+1
         t=t+dt
         tlmin=tlmin+dt
c         call coad1d (dt)
         call coad2d (dt)
         do i=1,nkt
            if (g(i).lt.0.) print *,i,g(i),'g(i).lt.0 !!!!!'
         enddo
         if (tlmin.ge.60.) then
            do ia=1,nka
            do jt=1,nkt
               if2(ia,jt)=f(ia,jt)
            enddo
            enddo
            tlmin=tlmin-60.
            lmin=lmin+1
            print *,'time in minutes:',lmin
            write (6) lday,lst,lmin,kf,ixm2,irh,izk,if2
            write (7,6110) (coldat(icol,i),i=1,10)
 6110 format (10f8.0)
c mass balance
            do i=1,nkt
               x1=max(x1,g(i))
               if (dabs(x1-g(i)).lt.1.d-9) imax=i
            enddo
            x0=0.
            x2=0.
            do ia=1,nka
            do jt=1,nkt
               x0=x0+f(ia,jt)*e(jt)
               x2=x2+f(ia,jt)*en(ia)*1.d09
            enddo
            enddo
            write (*,6000) x0*1000.,x2,imax
            write (*,6010) (x0/wini-1.)*100.,(x2/aerm-1.)*100.
 6000 format (1x,'water and aerosol mass:',2f16.8,i8)
 6010 format (1x,'mass loss in %:',8x,2f16.8)
            do i=2,nkt-1
               if (g(i).gt.g(i+1).and.g(i).gt.g(i-1)) then
                  if (g(i).gt.1.e-6) write (*,6005) g(i),i
               endif
            enddo
 6005 format (1x,'local maximum of water mass:',f16.8,i8)
         endif
      enddo
      return
      end

      subroutine coad2d (dt)
      implicit double precision (a-h,o-z)
      parameter (nf=100,n=nf+50,nkt=120,nka=120,nb=20,nka2=72)
      common /cb29/ c(nkt,nkt),elnr(nkt),ima(nkt,nkt),jmin
      common /cb50/ enw(nka),ew(nkt),rn(nka),rw(nka,nkt),en(nka),
     &              e(nkt),dew(nkt),rq(nka,nkt),kg,kgp,kd,kdp
      common /cb79/ xm1pl,dlnr,rl(nkt),ax,axn,i4
      common/kern/ck(nkt,nkt)
      common/courn/cn(nka,nka),icn(nka,nka)
      common /cb52/ fpi(nka),fpj(nka),f(nka,nkt),g(nkt),fs(nkt)
      dimension fk(nka),fkp(nka),fz(nka,nkt)
      data gmin,fmin /1.d-38,1.d-20/
c one-dimensional particle and mass distribution
      do jt=1,nkt
         fs(jt)=0.
         do ia=1,nka
            fz(ia,jt)=f(ia,jt)
            fs(jt)=fs(jt)+f(ia,jt)
         enddo
         g(jt)=fs(jt)*elnr(jt)
      enddo
c jtl, jtr: lower and upper boundary of water grid
      do jt=i4,nkt-1
         jtl=jt
         if (g(jt).gt.gmin) go to 2000
      enddo
 2000 continue
      do jt=nkt-1,jtl,-1
         jtr=jt
         if (g(jt).gt.gmin) go to 2010
      enddo
 2010 continue
      jtr=nkt-1
c initial total aerosol mass
      en0=0.
      do ia=1,nka
         x0=0.
         do jt=jtl,jtr
            x0=x0+f(ia,jt)
         enddo
         en0=en0+x0*en(ia)
      enddo
c i: loop over small particle distribution, i-drop collected by j-drop
      do i=jtl,jtr
c iali, iari: lower and upper boundary of aerosol grid in bin i
         do ia=1,nka
            iali=ia
            if (f(ia,i).gt.fmin) go to 2020
         enddo
         go to 2030
 2020    continue
         do ia=nka,iali,-1
            iari=ia
            if (f(ia,i).gt.fmin) go to 2040
         enddo
 2040    continue
c fsi: total particles in bin i 
         fsi=0.
         do ia=iali,iari
            fsi=fsi+f(ia,i)
         enddo
c fpi: fraction of particles at gridpoint (ia,i)
         do ia=iali,iari
            fpi(ia)=f(ia,i)/fsi
         enddo
c j: loop over large particle distribution, j-drop collects i-drop
         do j=i+1,jtr
c ialj, iarj: lower and upper boundary of aerosol grid in bin j
            do ia=1,nka
               ialj=ia
               if (f(ia,j).gt.fmin) go to 2050
            enddo
            go to 2030
 2050       continue
            do ia=nka,ialj,-1
               iarj=ia
               if (f(ia,j).gt.fmin) go to 2060
            enddo
 2060       continue
c fsj: total particles in bin j
            fsj=0.
            do ia=ialj,iarj
               fsj=fsj+f(ia,j)
            enddo
c fpj: fraction of particles at gridpoint (ia,j)
            do ia=ialj,iarj
               fpj(ia)=f(ia,j)/fsj
            enddo
c one-dimensional coalescence
            k=min(ima(i,j),nkt-1)
            kp=k+1
            x0=ck(i,j)*g(i)*g(j)*dt
            x0=min(x0,g(i)*e(j))
            if (j.ne.k) x0=min(x0,g(j)*e(i))
            gsi=x0/e(j)
            gsj=x0/e(i)
            gsk=gsi+gsj
            gjold=g(j)
            if (gsk.lt.gmin) go to 2030
            g(i)=g(i)-gsi
            g(j)=g(j)-gsj
            gk=g(k)+gsk
            if (gk.gt.gmin) then
               x1=dlog(g(kp)/gk+1.d-60)
               flux=gsk/x1*(dexp(0.5*x1)-dexp(x1*(0.5-c(i,j))))
               flux=min(flux,gk,gsk)
            else
               flux=0.
            endif
            fsk=(min(0.d0,g(k))+gsk-flux)/elnr(k)
            fskp=flux/elnr(kp)
            g(k)=gk-flux
            g(kp)=g(kp)+flux
c two-dimensional redistribution after coalescence
c counting of collision processes
c fac: number of i drops collected by one j drop
c            fac=max(1.d0,gsj/gjold)
c            call count (fac)
c update of particles in bins i and j 
            fsi=g(i)/elnr(i)
            do ia=iali,iari
               f(ia,i)=fsi*fpi(ia)
            enddo
            fsj=max(0.d0,(gjold-gsj))/elnr(j)
            do ia=ialj,iarj
               f(ia,j)=fsj*fpj(ia)
            enddo
c fk, fkp: new particles in bins k and k+1
            ialk=icn(iali,ialj)
            iark=min(nka,icn(iari,iarj)+1)
            do ia=ialk,iark
               fk(ia)=0.
               fkp(ia)=0.
            enddo
            do iai=iali,iari
            do iaj=ialj,iarj
               x1=fpi(iai)*fpj(iaj)
               iak=icn(iai,iaj)
               iakp=min(nka,iak+1)
               x2=cn(iai,iaj)
               x3=1.-x2
               fk(iak)=fk(iak)+fsk*x1*x3
               fk(iakp)=fk(iakp)+fsk*x1*x2
               fkp(iak)=fkp(iak)+fskp*x1*x3
               fkp(iakp)=fkp(iakp)+fskp*x1*x2
            enddo
            enddo
c update of two-dimensional particle distribution in k and k+1 bins
            do ia=ialk,iark
               f(ia,k)=f(ia,k)+fk(ia)
               f(ia,kp)=f(ia,kp)+fkp(ia)
            enddo
         enddo
 2030    continue
      enddo
c obtain aerosol mass balance by moving the particle distribution 
c in k and k+1 bins along the aerosol grid with upstream advection 
      en1=0.
      do ia=1,nka
         x0=0.
         do jt=jtl,jtr
            x1=f(ia,jt)
            x0=x0+x1
            if (x1.gt.fz(ia,jt)) then
               fz(ia,jt)=x1
               f(ia,jt)=0.
            else
               fz(ia,jt)=0.
            endif
         enddo
         en1=en1+x0*en(ia)
      enddo
      if (en0.gt.en1) then
         x0=0.
         do jt=jtl,jtr
         do ia=1,nka-1
            x0=x0+fz(ia,jt)*(en(ia+1)-en(ia))
         enddo
         enddo
         c0=min(1.d0,(en0-en1)/x0)
         x1=1.-c0
         do jt=jtl,jtr
            fz(nka,jt)=fz(nka,jt)+fz(nka-1,jt)*c0
            do ia=nka-1,2,-1
               fz(ia,jt)=fz(ia,jt)*x1+fz(ia-1,jt)*c0
            enddo
            fz(1,jt)=fz(1,jt)*x1
         enddo
      else
         x0=0.
         do jt=jtl,jtr
         do ia=2,nka
            x0=x0+fz(ia,jt)*(en(ia-1)-en(ia))
         enddo
         enddo
         c0=min(1.d0,(en0-en1)/x0)
         x1=1.-c0
         do jt=jtl,jtr
            fz(1,jt)=fz(1,jt)+fz(2,jt)*c0
            do ia=2,nka-1
               fz(ia,jt)=fz(ia,jt)*x1+fz(ia+1,jt)*c0
            enddo
            fz(nka,jt)=fz(nka,jt)*x1
         enddo
      endif
      do ia=1,nka
      do jt=jtl,jtr
         f(ia,jt)=f(ia,jt)+fz(ia,jt)
      enddo
      enddo
      return
      end

      subroutine coad1d (dt)
      implicit double precision (a-h,o-z)
      parameter (nkt=120,nka=120,nb=20,n=150)
      common /cb29/ c(nkt,nkt),elnr(nkt),ima(nkt,nkt),jmin
      common /cb50/ enw(nka),ew(nkt),rn(nka),rw(nka,nkt),en(nka),
     &              e(nkt),dew(nkt),rq(nka,nkt),kg,kgp,kd,kdp
      common /cb52/ fpi(nka),fpj(nka),f(nka,nkt),g(nkt),fs(nkt)
      common /cb79/ xm1pl,dlnr,rl(nkt),ax,axn,i4
      common/kern/ck(nkt,nkt)
      data gmin /1.d-38/
c g: one-dimensional mass distribution
      do jt=1,nkt
         g(jt)=f(1,jt)*elnr(jt)
      enddo
c lower and upper integration limit i0,i1
      do jt=i4,nkt-1
         jtl=jt
         if (g(jt).gt.gmin) go to 2010
      enddo
 2010 continue
      do jt=nkt-1,jtl,-1
         jtr=jt
         if (g(jt).gt.gmin) go to 2020
      enddo
 2020 continue
c one-dimensional coalescence
      do i=jtl,jtr
      do j=i+1,jtr
         k=ima(i,j)
         x0=ck(i,j)*g(i)*g(j)*dt
         x0=min(x0,g(i)*e(j))
         if (j.ne.k) x0=min(x0,g(j)*e(i))
         gs=x0*(1./e(i)+1./e(j))
         if (gs.ge.gmin) then
            g(i)=g(i)-x0/e(j)
            g(j)=g(j)-x0/e(i)
            g(k)=g(k)+gs
            x1=dlog((g(k+1)+1.d-60)/g(k))
            flux=gs/x1*(dexp(0.5*x1)-dexp(x1*(0.5d0-c(i,j))))
            flux=min(flux,g(k))        
            g(k)=g(k)-flux
            g(k+1)=g(k+1)+flux
         endif
      enddo
      enddo
      do jt=jtl,nkt
         f(1,jt)=g(jt)/elnr(jt)
      enddo
      return
      end

      subroutine count (fac)
      implicit double precision (a-h,o-z)
      common /datco/ coldat(70,10),icol
      if (fac.le.1.) coldat(icol,1)=coldat(icol,1)+1.
      if (fac.gt.1..and.fac.le.2.) coldat(icol,2)=coldat(icol,2)+1.
      if (fac.gt.2..and.fac.le.5.) coldat(icol,3)=coldat(icol,3)+1.
      if (fac.gt.5..and.fac.le.10.) coldat(icol,4)=coldat(icol,4)+1.
      if (fac.gt.10..and.fac.le.20.) coldat(icol,5)=coldat(icol,5)+1.
      if (fac.gt.20..and.fac.le.50.) coldat(icol,6)=coldat(icol,6)+1.
      if (fac.gt.50..and.fac.le.100.) coldat(icol,7)=coldat(icol,7)+1.
      if (fac.gt.100..and.fac.le.200.) coldat(icol,8)=coldat(icol,8)+1.
      if (fac.gt.200..and.fac.le.500.) coldat(icol,9)=coldat(icol,9)+1.
      coldat(icol,10)=coldat(icol,10)+1.
      return
      end

      subroutine fallg
      parameter(n=120)
      implicit double precision (a-h,o-z)
      dimension b(7),c(6)
      dimension rat(20),r0(15),ecoll(15,20)
      common/git/g(n),rdm(n),edm(n)
      common/veloc/winf(n)
      common/kern/ck(n,n)
      common/radcm/rr(n)
      data b/-0.318657e1,0.992696,-0.153193e-2,-0.987059e-3,
     &        -0.578878e-3,0.855176e-4,-0.327815e-5/
      data c/-0.500015e1,0.523778e1,-0.204914e1,0.475294,-0.542819e-1,
     &        0.238449e-2/
      data pi/3.141592654/
      eta=1.818e-4
      xlamb=6.62e-6
      rhow=1.
      rhoa=1.225e-3
      grav=980.665
      cunh=1.257*xlamb
      t0=273.15
      sigma=76.1-0.155*(293.15-t0)
      stok=2.*grav*(rhow-rhoa)/(9.*eta)
      stb=32.*rhoa*(rhow-rhoa)*grav/(3.*eta*eta)
      phy=sigma*sigma*sigma*rhoa*rhoa/((eta**4)*grav*(rhow-rhoa))
      py=phy**(1./6.)
c radius in [cm]-units
      do j=1,n
         rr(j)=rdm(j)*1.e-4
      enddo
      do j=1,n
         if (rr(j).le.1.e-3) then
            winf(j)=stok*(rr(j)*rr(j)+cunh*rr(j))
         elseif (rr(j).gt.1.e-3.and.rr(j).le.5.35e-2) then
            x=log(stb*rr(j)*rr(j)*rr(j))
            y=0.
            do i=1,7
               y=y+b(i)*(x**(i-1))
            enddo
            xrey=(1.+cunh/rr(j))*exp(y)
            winf(j)=xrey*eta/(2.*rhoa*rr(j))
         elseif (rr(j).gt.5.35e-2) then
            bond=grav*(rhow-rhoa)*rr(j)*rr(j)/sigma
            if (rr(j).gt.0.35) bond=grav*(rhow-rhoa)*0.35*0.35/sigma
            x=log(16.*bond*py/3.)
            y=0.
            do i=1,6
               y=y+c(i)*(x**(i-1))
            enddo
            xrey=py*exp(y)
            winf(j)=xrey*eta/(2.*rhoa*rr(j))
            if (rr(j).gt.0.35)  winf(j)=xrey*eta/(2.*rhoa*0.35)
         endif
      enddo
      return
      end

      subroutine trkern(isw)
      implicit double precision (a-h,o-z)
      parameter (nkt=120,nka=120,nb=20,n=nkt,nktt=80,nkaa=40)
      common /cb29/ c(nkt,nkt),elnr(nkt),ima(nkt,nkt),jmin
      common /cb50/ enw(nka),ew(nkt),rn(nka),rw(nka,nkt),en(nka),
     &              edm(nkt),dew(nkt),rq(nka,nkt),kg,kgp,kd,kdp
      common /cb52/ fpi(nka),fpj(nka),f(nka,nkt),g(nkt),fs(nkt)
      common /cb79/ xm1pl,dlnr,rdm(nkt),ax,axn,i4
      common/courn/cn(nka,nka),icn(nka,nka)
      common/kern/ck(nkt,nkt)
      common/veloc/winf(nkt)
c      common/kern/ck(n,n),ec(n,n)
      common/radcm/rr(nkt)
      common /co2d/aerm,wini,fogtype,filename
      dimension cck(nkt,nkt)

      data pi/3.141592654/
      if (isw.eq.0) then
c LONG-KERNEL
         do j=1,n
         do i=1,j
            if(rdm(j).le.50.) then
               effi=4.5d-4*rdm(j)*rdm(j)*
     &              (1.d0-3.d0/(max(3.d0,dble(rdm(i)))+1.d-2))
            else
               effi=1.d0
            endif
            cck(j,i)=pi*(rr(j)+rr(i))*(rr(j)+rr(i))*effi*
     &              abs(winf(j)-winf(i))
            cck(i,j)=cck(j,i)
          enddo
          enddo
      elseif (isw.eq.1) then
c REAL-KERNEL
c         call effic
c         do j=1,n
c         do i=1,j
c            cck(j,i)=pi*(rr(j)+rr(i))*(rr(j)+rr(i))*ec(j,i)*
c     &              abs(winf(j)-winf(i))
c            cck(i,j)=cck(j,i)
c         enddo
c         enddo
      else
c GOLOVIN-KERNEL
         do j=1,n
         do i=1,j
            cck(j,i)=1.5*(edm(j)+edm(i))
            cck(i,j)=cck(j,i)
         enddo
         enddo
      endif
      do i=1,n
      do j=1,n
         jm=max0(j-1,1)
         im=max0(i-1,1)
         jp=min0(j+1,n)
         ip=min0(i+1,n)
         ck(i,j)=0.125*(cck(i,jm)+cck(im,j)+cck(ip,j)+cck(i,jp))
     &          +.5*cck(i,j)
         if (i.eq.j) ck(i,j)=0.5*ck(i,j)
      enddo
      enddo
      return
      end


      subroutine courant (fogtype)
      include 'parameter.h'
      include 'cb29.h'
      include 'cb50.h'
      include 'cb52.h'
      include 'cb79.h'
      dimension ckk(nkt,nkt)
      character *1 fogtype
      character *7 fname
c courant numbers of water grid
      do i=1,nkt
         do j=1,nkt
            x0=e(i)+e(j)
            if (x0.ge.e(nkt)-1.d-10) then
               c(i,j)=0.d0
               ima(i,j)=nkt
            else
               do k=j,nkt
                  k0=k-1
                  if (e(k).ge.x0.and.e(k0).lt.x0) then
                     c(i,j)=dlog(x0/e(k0))/(3.d0*dlnr)
                     if (c(i,j).gt.1.-1.d-08) then
                        c(i,j)=0.
                        k0=k
                     endif
                     ima(i,j)=k0
                     goto 2000
                  endif
               enddo
            endif
 2000    continue
         enddo
      enddo
c elnr: conversion from particle to mass distribution
c jmin: minimum bin of collected particles
      do jt=1,nkt
         elnr(jt)=e(jt)/dlnr
         if (rl(jt).lt.1.) jmin=jt
      enddo
c courant numbers of aerosol grid
      do iai=1,nka
         do iaj=1,nka
            x0=en(iai)+en(iaj)
            if (x0.gt.en(nka)) then
               cn(iai,iaj)=0.d0
               icn(iai,iaj)=nka
            else
               do iak=iai,nka-1
                  if (x0.gt.en(iak).and.x0.le.en(iak+1)) then
c                  cn(iai,iaj)=(x0-en(iak))/(en(iak+1)-en(iak))
c                  cn(iai,iaj)=min(cn(iai,iaj),1.d0)
                     cn(iai,iaj)=dlog(x0/en(iak))
     &                           /dlog(en(iak+1)/en(iak))
                     icn(iai,iaj)=iak
                     go to 2010
                  endif
               enddo
 2010          continue
            endif
         enddo
      enddo
c output for offline calculation of coalescence
c      do ia=1,nkt
c      do jt=1,nkt
cccccc         ckk(ia,jt)=ck(2,ia,jt)
c         ckk(ia,jt)=ck(ia,jt)
c      enddo
c      enddo
c      fname='ci .out'
c      fname(3:3)=fogtype
c      open (77, file=fname,status='new',form='unformatted')
c      write (77) cn,ckk,c,elnr,rn,en,e,rl,dlnr,ax,axn,ima,jmin,icn
c      close (77)
      return
      end

