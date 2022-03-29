      subroutine gibbs(iseed)
      implicit real*8 (a-h, o-z)
	integer, parameter :: ni = 11394, nj = 50
	integer, parameter :: nc = 5, nb = 3, nl = 2
			
	real*8 yi(ni), zi(ni,nc), xij(ni,nj,nb)    
      real*8 wi(ni), taui(ni), vnu 
	real*8 betaa(nc), alphaa(nb,nl)
      real*8 azeta(nl), agamma
      real*8 gj(nj), weight(nl), theta(nl)
      
      real*8 sigma_betaa, sigma_theta, a0, b0

      real*8 prob(nj,nl)
      
      common /vyi/yi
      common /vzi/zi
      common /vxij/xij
      
      common /vwi/wi
      common /vtaui/taui
      common /vvnu/vnu

      common /vbetaa/betaa
      
      common /valphaa/alphaa
      common /vazeta/azeta
      common /vagamma/agamma
      
      common /vgj/gj
      common /vweight/weight
      common /vtheta/theta
            
      common /vsigma_betaa/sigma_betaa
      common /vsigma_theta/sigma_theta
      common /va0/a0
      common /vb0/b0
      
      common /vprob/prob
      
      call gibbs_gj(iseed) 
      call gibbs_theta(iseed)          
      
      call gibbs_wi(iseed)          
      call gibbs_taui(iseed)          
      
      call gibbs_alphaa(iseed)                
      call gibbs_betaa(iseed)          
                  
      call gibbs_agamma(iseed)               
      call gibbs_azeta(iseed)         
      
      end subroutine
          
	subroutine gibbs_gj(iseed)  
      implicit real*8 (a-h, o-z)
	integer, parameter :: ni = 11394, nj = 50
	integer, parameter :: nc = 5, nb = 3, nl = 2
		
	real*8 yi(ni), zi(ni,nc), xij(ni,nj,nb)    
      real*8 vnu, betaa(nc), alphaa(nb,nl)
      real*8 gj(nj), weight(nl)
     
      real*8 prob(nj,nl)
      
      real*8 zb, xa, fxij 
      real*8 cdf, cprob(nl)
      real*8 amax, u, summ
      
      common /vyi/yi
      common /vzi/zi
      common /vxij/xij
      
      common /vvnu/vnu

      common /vbetaa/betaa
      common /valphaa/alphaa
      
      common /vgj/gj
      common /vweight/weight
      
      common /vprob/prob
      
      external drnunf, dtdf
                  
      do j = 1, nj
          
          do l = 1, nl
              prob(j,l) = dlog(weight(l))
          enddo
          do i = 1, ni
            
              zb = 0.d0
              do jj = 1, nc
                  zb = zb + zi(i,jj)*betaa(jj)
              enddo
                                                           
              fxij = 0.d0               
              do j1 = 1, nj
                  if (j1 .ne. j) then
                  if (xij(i,j1,1) .ne. 0.d0) then
              
                      l = nint(gj(j1))
                                    
                      xa = 0.d0
                      do j2 = 1, nb
                          xa = xa + xij(i,j1,j2)*alphaa(j2,l)
                      enddo
                      
                      fxij = fxij + xa
                      
                  endif
                  endif
              enddo
          
              do l = 1, nl
                            
                  xa = 0.d0
                  do jj = 1, nb
                      xa = xa + xij(i,j,jj)*alphaa(jj,l)
                  enddo
                      
                  cdf = dtdf(zb + fxij + xa, vnu)
                  
                  if (yi(i) .eq. 1.d0) then
                                            
                      prob(j,l) = prob(j,l) + dlog(cdf)
                      
                  else

                      prob(j,l) = prob(j,l) + dlog(1.d0 - cdf)
                      
                  endif
              
              enddo
                        
          enddo
         	
          amax = prob(j,1)
          do l = 2, nl
              if (amax .lt. prob(j,l)) then
                  amax = prob(j,l)
              endif
          enddo       
                
          summ = 0.d0
          do l = 1, nl
              prob(j,l) = dexp(prob(j,l) - amax)
              summ = summ + prob(j,l)
          enddo
          do l = 1, nl
              prob(j,l) = prob(j,l)/summ
          enddo
        
          cprob(1) = prob(j,1)
          do l = 2, nl
              cprob(l) = cprob(l-1) + prob(j,l)
          enddo

          call rnset(iseed)
          u = drnunf()
          call rnget(iseed)
                		                        
          if (u .le. cprob(1)) then
              gj(j) = 1.d0
          endif
          do l = 2, nl
              if ((u .gt. cprob(l-1)) 
     +             .and. (u .le. cprob(l))) then
                  gj(j) = dfloat(l)
              endif
          enddo
          
      enddo
                  
      end subroutine	   
                   

      subroutine gibbs_theta(iseed) 
      implicit real*8 (a-h, o-z)
	integer, parameter :: ni = 11394, nj = 50
	integer, parameter :: nc = 5, nb = 3, nl = 2

      real*8 gj(nj), weight(nl), theta(nl)
      real*8 sigma_theta
      
	real*8 start(1), xmin(1), ynewlo, reqmin, step(1)
	integer konvge, kcount, icount, numres, ifault, nopt

	real*8 bold, anew, amean, asigma
	real*8 rv, u, bpdf, apdf, ratio
      
      integer ldum
      real*8 ttheta(nl), sumw
      real*8 suma, sumb, sumc, sumd, summ
            
      common /vgj/gj
      common /vweight/weight
      common /vtheta/theta
                  
      common /vsigma_theta/sigma_theta
      
	common /vldum/ldum
	common /vttheta/ttheta
      
      external ftheta, drnnof, drnunf
            
      ttheta(1) = 0.d0      
      do l = 2, nl
          if (l .eq. 2) then
              ttheta(l) = dlog(-theta(l))
          else
              ttheta(l) = dlog(theta(l-1) - theta(l))
          endif
      enddo
      
      do l = 2, nl
          
          ldum = l
          
          bold = ttheta(l)
                                  
          nopt = 1
	    reqmin = 1.0d-10 ; konvge = 5 ; kcouns = 1000
	    step(1) = 0.2d0 ; start(1) = bold
          call nelmin(ftheta, nopt, start, xmin, ynewlo,  
     +                reqmin, step, konvge, kcount, 
     +                icount, numres, ifault)
	    amean = xmin(1)
          ttheta(l) = amean
                             
          theta(1) = 0.d0      
          do ll = 2, nl
              if (ll .eq. 2) then
                  theta(ll) = -dexp(ttheta(ll))
              else
                  theta(ll) = theta(ll-1) - dexp(ttheta(ll))
              endif
          enddo
          
          summ = 0.d0
          do j = 1, nj          

              suma = 1.d0
              do ll = 2, nl
                  suma = suma + dexp(theta(ll))                  
              enddo
              
              sumb = 0.d0; sumc = 0.d0 
              do ll = l, nl
                  
                  if (nint(gj(j)) .eq. ll) then
                      sumb = sumb + 1.d0
                  endif
                      
                  sumc = sumc + dexp(theta(ll))
                                    
              enddo

              summ = summ 
     +             - sumb*dexp(ttheta(l))
     +             + ( dexp(ttheta(l))
     +                 *(1.d0 - dexp(ttheta(l)))
     +                 *sumc*suma 
     +                 + (dexp(ttheta(l))*sumc)**2 )
     +               /suma**2
                            
          enddo
          
          sumd = 0.d0 
          do ll = l, nl
                                    
              sumd = sumd
     +             + dexp(ttheta(l))*theta(ll)/sigma_theta
     +             - dexp(2.d0*ttheta(l))/sigma_theta
                  
          enddo  
          
          asigma = summ + sumd
          
          asigma = -1.d0/asigma*1.5d0
          if (asigma .lt. 0.0d0) asigma = -asigma      
                 
          bpdf = -ftheta(bold) 
     +         + (bold - amean)**2/(2.d0*asigma)
        
          do ii = 1, 100
	  
              call rnset(iseed)
              rv = drnnof()
              call rnget(iseed)
              anew = amean + rv*dsqrt(asigma)

              apdf = -ftheta(anew) 
     +             + (anew - amean)**2/(2.d0*asigma)
              ratio = apdf - bpdf 
     
              if (ratio .ge. 0.0d0) then
                  bpdf = apdf
                  bold = anew
              else
                  call rnset(iseed)
                  u = drnunf()
                  call rnget(iseed)
                  if (dlog(u) .le. ratio) then
                      bpdf = apdf
                      bold = anew
                  endif
              endif
		
          enddo     
                    
          ttheta(l) = bold
                                   
      enddo
      
      theta(1) = 0.d0      
      do l = 2, nl
          if (l .eq. 2) then
              theta(l) = -dexp(ttheta(l))
          else
              theta(l) = theta(l-1) - dexp(ttheta(l))
          endif
      enddo
                    
      sumw = 0.d0
      do l = 1, nl
          sumw = sumw + dexp(theta(l))
      enddo
      do l = 1, nl
          weight(l) = dexp(theta(l))/sumw
      enddo
      
      end subroutine                   
             
      real*8 function ftheta(star)
      implicit real*8 (a-h, o-z)
	integer, parameter :: ni = 11394, nj = 50
	integer, parameter :: nc = 5, nb = 3, nl = 2

      real*8 gj(nj)
      real*8 weight(nl), theta(nl)
      real*8 sigma_theta

      integer ldum
      real*8 ttheta(nl)

      real*8 star, sumw, sum1, sum2, pdf
      
      common /vgj/gj
                  
      common /vsigma_theta/sigma_theta
      
	common /vldum/ldum
	common /vttheta/ttheta
            
      l = ldum
                   
      ttheta(l) = star
                   
      theta(1) = 0.d0      
      do ll = 2, nl
          if (ll .eq. 2) then
              theta(ll) = -dexp(ttheta(ll))
          else
              theta(ll) = theta(ll-1) - dexp(ttheta(ll))
          endif
      enddo
                
      sumw = 0.d0
      do ll = 1, nl
          sumw = sumw + dexp(theta(ll))
      enddo
      do ll = 1, nl
          weight(ll) = dexp(theta(ll))/sumw
      enddo
                             
      sum1 = 0.d0
      do j = 1, nj          
          ll = nint(gj(j))
          sum1 = sum1 + dlog(weight(ll))
      enddo
            
      sum2 = 0.d0
      do ll = 2, nl
          sum2 = sum2 
     +         + ttheta(ll) 
     +         - theta(ll)**2/(2.d0*sigma_theta)
      enddo
      
      pdf = sum1 + sum2 
                          
      ftheta = -pdf
                        
      end function 
         
      
      subroutine gibbs_agamma(iseed)
      implicit real*8 (a-h, o-z)
	integer, parameter :: ni = 11394, nj = 50
	integer, parameter :: nc = 5, nb = 3, nl = 2

	real*8 alphaa(nb,nl), agamma
      real*8 a0, b0
      
	real*8 start(1), xmin(1), ynewlo, reqmin, step(1)
	integer konvge, kcount, icount, numres, ifault, nopt

	real*8 bold, anew, amean, asigma
	real*8 rv, u, bpdf, apdf, ratio
      
      real*8 summ, temp
      
      common /valphaa/alphaa
      common /vagamma/agamma
                        
      common /va0/a0
      common /vb0/b0
            
      external fagamma, drnnof, drnunf
       
      bold = dlog(agamma)

      nopt = 1
      reqmin = 1.0d-10 ; konvge = 5 ; kcouns = 1000
      step(1) = 0.2d0 ; start(1) = bold
      call nelmin(fagamma, nopt, start, xmin, ynewlo,  
     +            reqmin, step, konvge, kcount, 
     +            icount, numres, ifault)
      amean = xmin(1)
      agamma = dexp(amean)
      
      summ = 0.d0
      do l = 1, nl                              
          
          temp = 0.d0
          do jj = 1, nb
              temp = temp + alphaa(jj,l)**2
          enddo          
          
          summ = summ 
     +         - dsqrt(temp)*dexp(amean/2.d0)/4.d0
          
      enddo
            
      asigma = summ - b0*dexp(amean)
                
      asigma = -1.d0/asigma*1.5d0       
      if (asigma .lt. 0.0d0) asigma = -asigma      
                 
      bpdf = -fagamma(bold) 
     +     + (bold - amean)**2/(2.d0*asigma)
        
      do ii = 1, 100
	  
          call rnset(iseed)
          rv = drnnof()
          call rnget(iseed)
          anew = amean + rv*dsqrt(asigma)

          apdf = -fagamma(anew) 
     +         + (anew - amean)**2/(2.d0*asigma)
          ratio = apdf - bpdf 
     
          if (ratio .ge. 0.0d0) then
              bpdf = apdf
              bold = anew
          else
              call rnset(iseed)
              u = drnunf()
              call rnget(iseed)
              if (dlog(u) .le. ratio) then
                  bpdf = apdf
                  bold = anew
              endif
          endif
		
      enddo     

      agamma = dexp(bold)      
          
      end subroutine                   
   
      real*8 function fagamma(star)
      implicit real*8 (a-h, o-z)
	integer, parameter :: ni = 11394, nj = 50
	integer, parameter :: nc = 5, nb = 3, nl = 2

	real*8 alphaa(nb,nl)
      real*8 a0, b0
      
      real*8 star, summ, temp, pdf
      
      common /valphaa/alphaa
                        
      common /va0/a0
      common /vb0/b0
                  
      summ = 0.d0
      do l = 1, nl                              
          
          temp = 0.d0
          do jj = 1, nb
              temp = temp + alphaa(jj,l)**2
          enddo          
          
          summ = summ 
     +         + dfloat(nb)*star/2.d0
     +         - dsqrt(temp)*dexp(star/2.d0)
          
      enddo
      
      pdf = summ + a0*star - b0*dexp(star)

      fagamma = -pdf
    
      end function

      
      subroutine gibbs_azeta(iseed)
      implicit real*8 (a-h, o-z)
	integer, parameter :: ni = 11394, nj = 50
	integer, parameter :: nc = 5, nb = 3, nl = 2

	real*8 alphaa(nb,nl)
      real*8 azeta(nl), agamma 
      
      real*8 pp, aa, bb, temp, rv
      
      common /valphaa/alphaa
      common /vazeta/azeta
      common /vagamma/agamma
                                    
      external rgenGIG
         
      do l = 1, nl
                    
          temp = 0.d0
          do jj = 1, nb
              temp = temp + alphaa(jj,l)**2
          enddo
          
          pp = -1.d0/2.d0
          aa = agamma*temp
          bb = 1.d0
          
          rv = rgenGIG(pp,aa,bb,iseed)
          azeta(l) = 1.d0/rv
          
      enddo
                    
      end subroutine                   
        
   
      subroutine gibbs_alphaa(iseed)
      implicit real*8 (a-h, o-z)
	integer, parameter :: ni = 11394, nj = 50
	integer, parameter :: nc = 5, nb = 3, nl = 2

	real*8 zi(ni,nc), xij(ni,nj,nb)    
      real*8 wi(ni), taui(ni)
	real*8 alphaa(nb,nl) 
      real*8 azeta(nl), agamma
      real*8 gj(nj)
                
      real*8 sigma_betaa
      
      real*8 xstar(nb), wstar, xa, fxij
      real*8 sumxx(nb,nb), sumxw(nb)
      real*8 sumzz(nc,nc), sumzx(nc,nb), sumzw(nc)
      real*8 AA(nc,nc), AAi(nc,nc)
      real*8 BB(nb,nc), CC(nb,nb), DD(nb)
      
      real*8 hmean(nb), hsigma(nb,nb)
      real*8 amean(nb), asigma(nb,nb)
	real*8 tol, rsig(nb,nb), rv(nb)
      
      common /vzi/zi
      common /vxij/xij
      
      common /vwi/wi
      common /vtaui/taui

      common /valphaa/alphaa
      common /vazeta/azeta
      common /vagamma/agamma
      
      common /vgj/gj
            
      common /vsigma_betaa/sigma_betaa
      
      external dlinrg, dblinf, dmach, dchfac, drnmvn
            
      do l = 1, nl
                                                      
          do j1 = 1, nb
              
              sumxw(j1) = 0.d0
              
              do j2 = 1, nb
                  sumxx(j1,j2) = 0.d0
              enddo
                      
          enddo
              
          do j1 = 1, nc      
        
              sumzw(j1) = 0.d0
              
              do j2 = 1, nc           
                  sumzz(j1,j2) = 0.d0
              enddo

              do j2 = 1, nb           
                  sumzx(j1,j2) = 0.d0
              enddo
                  
          enddo             
          
          do i = 1, ni
                                  
              do jj = 1, nb
                  xstar(jj) = 0.d0
              enddo              
              fxij = 0.d0               
              do j = 1, nj              
                  if (xij(i,j,1) .ne. 0.d0) then
                  
                      ll = nint(gj(j))
                
                      if (ll .eq. l) then
                      
                          do jj = 1, nb
                              xstar(jj) = xstar(jj) + xij(i,j,jj)
                          enddo
                      
                      else                          
                          
                          xa = 0.d0
                          do jj = 1, nb
                              xa = xa + xij(i,j,jj)*alphaa(jj,ll)
                          enddo
                      
                          fxij = fxij + xa
                          
                      endif  
                      
                  endif  
              enddo
                            
              wstar = wi(i) - fxij
                                                    
              do j1 = 1, nb
                      
                  sumxw(j1) = sumxw(j1) + taui(i)*xstar(j1)*wstar
                      
                  do j2 = 1, nb
                  
                      sumxx(j1,j2) = sumxx(j1,j2) 
     +                             + taui(i)*xstar(j1)*xstar(j2)
                      
                  enddo
                      
              enddo
              
              do j1 = 1, nc      
        
                  sumzw(j1) = sumzw(j1) + taui(i)*zi(i,j1)*wstar
            
                  do j2 = 1, nc           
                
                      sumzz(j1,j2) = sumzz(j1,j2) 
     +                             + taui(i)*zi(i,j1)*zi(i,j2)

                  enddo

                  do j2 = 1, nb           
                
                      sumzx(j1,j2) = sumzx(j1,j2) 
     +                             + taui(i)*zi(i,j1)*xstar(j2)

                  enddo
                  
              enddo             
                                                  
          enddo
            
          do j1 = 1, nc                           
              do j2 = 1, nc           
                  AA(j1,j2) = sumzz(j1,j2)
              enddo                  
              AA(j1,j1) = AA(j1,j1) + 1.d0/sigma_betaa
          enddo      
                    
          call dlinrg(nc, AA, nc, AAi, nc)

          do j1 = 1, nb                           
              do j2 = 1, nc
                  BB(j1,j2) = 0.d0
                  do jj = 1, nc
                      BB(j1,j2) = BB(j1,j2) 
     +                          + sumzx(jj,j1)*AAi(jj,j2)
                  enddo
              enddo                  
          enddo             

          do j1 = 1, nb 
              
              do j2 = 1, nb
                  CC(j1,j2) = 0.d0
                  do jj = 1, nc
                      CC(j1,j2) = CC(j1,j2) 
     +                          + BB(j1,jj)*sumzx(jj,j2)
                  enddo
              enddo        
              
              DD(j1) = 0.d0
              do j2 = 1, nc
                  DD(j1) = DD(j1) + BB(j1,j2)*sumzw(j2)
              enddo
              
          enddo             
          
          do j1 = 1, nb     
        
              hmean(j1) = sumxw(j1) - DD(j1)
            
              do j2 = 1, nb
                
                  hsigma(j1,j2) = sumxx(j1,j2) 
     +                          -  CC(j1,j2)

              enddo
              
              hsigma(j1,j1) = hsigma(j1,j1) 
     +                      + agamma/azeta(l)
              
          enddo             
          
          call dlinrg(nb, hsigma, nb, asigma, nb)
	
          do j1 = 1, nb
              amean(j1) = 0.d0
              do j2 = 1, nb
                  amean(j1) = amean(j1) 
     +                      + asigma(j1,j2)*hmean(j2)
              enddo
          enddo

          tol = 100.d0*dmach(4)
          call dchfac(nb, asigma, nb, tol, irank, rsig, nb)
          call rnset(iseed)
          call drnmvn(1, nb, rsig, nb, rv, 1)
          call rnget(iseed)      
	
          do jj = 1, nb
              alphaa(jj,l)= amean(jj) + rv(jj)
          enddo  
          
      enddo
       
      end subroutine                   
                 

      subroutine gibbs_betaa(iseed)
      implicit real*8 (a-h, o-z)
	integer, parameter :: ni = 11394, nj = 50
	integer, parameter :: nc = 5, nb = 3, nl = 2
		
	real*8 zi(ni,nc), xij(ni,nj,nb)    
      real*8 wi(ni), taui(ni)
	real*8 betaa(nc), alphaa(nb,nl)
      real*8 gj(nj)
            
      real*8 sigma_betaa
      
      real*8 hmean(nc), hsigma(nc,nc)
      real*8 amean(nc), asigma(nc,nc)
	real*8 tol, rsig(nc,nc), rv(nc)
      
      real*8 wstar, xa, fxij
      
      common /vzi/zi
      common /vxij/xij
      
      common /vwi/wi
      common /vtaui/taui

      common /vbetaa/betaa
      common /valphaa/alphaa
      
      common /vgj/gj
            
      common /vsigma_betaa/sigma_betaa
      
      external dmach, dchfac, drnmvn
      
      do j1 = 1, nc
          
          hmean(j1) = 0.d0
          
          do j2 = 1, nc
              hsigma(j1,j2) = 0.d0
          enddo
          hsigma(j1,j1) = 1.d0/sigma_betaa
                        
      enddo
      
      do i = 1, ni
                
          fxij = 0.d0               
          do j = 1, nj  
              if (xij(i,j,1) .ne. 0.d0) then
              
                  l = nint(gj(j))
                                    
                  xa = 0.d0
                  do jj = 1, nb
                      xa = xa + xij(i,j,jj)*alphaa(jj,l)
                  enddo
                      
                  fxij = fxij + xa
                  
              endif
          enddo
          
          wstar = wi(i) - fxij
          
          do j1 = 1, nc      
        
              hmean(j1) = hmean(j1) + taui(i)*zi(i,j1)*wstar
            
              do j2 = 1, nc           
                
                  hsigma(j1,j2) = hsigma(j1,j2) 
     +                          + taui(i)*zi(i,j1)*zi(i,j2)

              enddo

          enddo             
                              		                        
      enddo  
                           
	call dlinrg(nc, hsigma, nc, asigma, nc)
	
	do j1 = 1, nc
	    amean(j1) = 0.d0
		do j2 = 1, nc
		    amean(j1) = amean(j1) + asigma(j1,j2)*hmean(j2)
		enddo
      enddo

	tol = 100.d0*dmach(4)
	call dchfac(nc, asigma, nc, tol, irank, rsig, nc)
	call rnset(iseed)
	call drnmvn(1, nc, rsig, nc, rv, 1)
	call rnget(iseed)   
		
	do jj = 1, nc
          betaa(jj) = amean(jj) + rv(jj)
      enddo      
      
      end subroutine            
   

      subroutine gibbs_taui(iseed)
      implicit real*8 (a-h, o-z)
	integer, parameter :: ni = 11394, nj = 50
	integer, parameter :: nc = 5, nb = 3, nl = 2
		
	real*8 zi(ni,nc), xij(ni,nj,nb)    
      real*8 wi(ni), taui(ni), vnu 
	real*8 betaa(nc), alphaa(nb,nl)
      real*8 gj(nj)
            
      real*8 wstar, zb, xa, fxij
      real*8 shape, scale, rv
      
      common /vzi/zi
      common /vxij/xij
      
      common /vwi/wi
      common /vtaui/taui
      common /vvnu/vnu

      common /vbetaa/betaa
      common /valphaa/alphaa
      
      common /vgj/gj
                  
      external drngam
      
      do i = 1, ni
      
          zb = 0.d0
          do jj = 1, nc
              zb = zb + zi(i,jj)*betaa(jj)
          enddo
          
          fxij = 0.d0               
          do j = 1, nj  
              if (xij(i,j,1) .ne. 0.d0) then
              
                  l = nint(gj(j))
                                    
                  xa = 0.d0
                  do jj = 1, nb
                      xa = xa + xij(i,j,jj)*alphaa(jj,l)
                  enddo
                      
                  fxij = fxij + xa
                  
              endif
          enddo
                                                          
          wstar = wi(i) - zb - fxij
                    
		shape = (vnu + 1.d0)/2.d0
		scale = (vnu + wstar**2)/2.d0
          
          call rnset(iseed)
          call drngam(1,shape,rv)
          call rnget(iseed)

          taui(i) = rv/scale      
                    		                        
      enddo  
                         
      end subroutine      
        
      
      subroutine gibbs_wi(iseed)
      implicit real*8 (a-h, o-z)
	integer, parameter :: ni = 11394, nj = 50
	integer, parameter :: nc = 5, nb = 3, nl = 2
		
	real*8 yi(ni), zi(ni,nc), xij(ni,nj,nb)    
      real*8 wi(ni), taui(ni)
	real*8 betaa(nc), alphaa(nb,nl)
      real*8 gj(nj)
      
      real*8 zb, xa, fxij
      real*8 wmean, wsigma, trim, rv       
      
      logical la, lb
            
      common /vyi/yi
      common /vzi/zi
      common /vxij/xij
      
      common /vwi/wi
      common /vtaui/taui

      common /vbetaa/betaa
      common /valphaa/alphaa
      
      common /vgj/gj
                  
      external ytuvn
      
      do i = 1, ni
      
          zb = 0.d0
          do jj = 1, nc
              zb = zb + zi(i,jj)*betaa(jj)
          enddo
          
          fxij = 0.d0               
          do j = 1, nj  
              if (xij(i,j,1) .ne. 0.d0) then
              
                  l = nint(gj(j))
                                    
                  xa = 0.d0
                  do jj = 1, nb
                      xa = xa + xij(i,j,jj)*alphaa(jj,l)
                  enddo
                      
                  fxij = fxij + xa
                  
              endif
          enddo
                                            
          wmean = zb + fxij
          wsigma = 1.d0/taui(i)

		trim = -wmean/dsqrt(wsigma)
          
		if (yi(i) .eq. 0.d0) then
		    la = .true. ; lb = .false.
          else
		    la = .false. ; lb = .true.
		endif

		call rnset(iseed)
		rv = ytuvn(trim, trim, la, lb, iseed)
		call rnget(iseed)
          
		wi(i) = wmean + rv*dsqrt(wsigma)                             
          	          
      enddo  
                
      end subroutine  