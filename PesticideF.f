	program main
      implicit real*8 (a-h, o-z)
	integer, parameter :: ni = 11394, nj = 50
	integer, parameter :: nc = 5, nb = 3, nl = 2
      integer, parameter :: iiniter = 1000, initer = 20000
      integer, parameter :: niter = 100000, niter1 = 20000
	      
	real*8 cancer, base(nc-1)
      real*8 durat(nj), ddurat
      real*8 zmean(nc), zstd(nc)
	
	real*8 yi(ni), zi(ni,nc), xij(ni,nj,nb)    
      real*8 wi(ni), taui(ni), vnu 
	real*8 betaa(nc), alphaa(nb,nl)
      real*8 azeta(nl), agamma
      real*8 gj(nj), weight(nl), theta(nl)
      
      real*8 sigma_betaa, sigma_theta, a0, b0
      
      integer ngj(nl)
      real*8 sumw

      real*8 cngj(nj,nl)
      real*8 prob(nj,nl)
                        
      real*8 delta(nj,nb)
	real*8 mean_betaa(nc)
      real*8 mean_delta(nj,nb)      
	real*8 bardic, dicbar, pd, dic      
      
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
            
      external gen_dic
            
    	open(unit = 5, file = 'nAHSData.txt') 
     	open(unit = 6, file = 'PesticideF_Initial.txt') 
     
    	open(unit = 11, file = 'PesticideF_Output1.txt') 
    	open(unit = 12, file = 'PesticideF_Output2.txt') 
    	open(unit = 13, file = 'PesticideF_Output3.txt') 
    	open(unit = 14, file = 'PesticideF_Output4.txt') 
    	open(unit = 15, file = 'PesticideF_Output5.txt') 
    	open(unit = 16, file = 'PesticideF_Output6.txt') 
    	open(unit = 17, file = 'PesticideF_Output7.txt') 
    	open(unit = 18, file = 'PesticideF_Output8.txt') 
      
      iseed = 9999999
                
      do jj = 1, nc
          zmean(jj) = 0.d0; zstd(jj) = 0.d0
      enddo   
	do ii = 1, ni
		             
	    read(5,*) i, cancer, base, durat
                    
	    yi(i) = cancer
          
          zi(i,1) = 1.d0
          do jj = 2, nc
              zi(i,jj) = base(jj-1)
              zmean(jj) = zmean(jj) + zi(i,jj)/dfloat(ni)
              zstd(jj) = zstd(jj) + zi(i,jj)**2
          enddo
          
          do j = 1, nj
              xij(i,j,1) = durat(j)
              xij(i,j,2) = durat(j)**2
              xij(i,j,3) = durat(j)**3
          enddo
                    
      enddo 
               
      do jj = 2, nc  
          temp = (zstd(jj) - dfloat(ni)*zmean(jj)**2)
     +            /dfloat(ni-1)
          zstd(jj) = dsqrt(temp)
      enddo

      do i = 1, ni
          do jj = 2, nc  
              zi(i,jj) = (zi(i,jj) - zmean(jj))/zstd(jj)
          enddo
      enddo
      
c     set hyper-parameters
                  
      sigma_betaa = 1000.d0
      sigma_theta = 1000.d0
      a0 = 1.0d0 ; b0 = 0.1d0
      vnu = 7.d0 
                       
c     set initial values    
                      
      read(6,*) betaa
      read(6,*) alphaa
      read(6,*) azeta
      read(6,*) agamma
      read(6,*) theta
      
      sumw = 0.d0
      do l = 1, nl
          sumw = sumw + dexp(theta(l))
      enddo
      do l = 1, nl
          weight(l) = dexp(theta(l))/sumw
      enddo
                 
      do i = 1, ni
          wi(i) = 0.d0
          taui(i) = 1.d0
      enddo                     
      
      do j = 1, nj
          gj(j) = 1.d0
      enddo

c     Run Gibbs
      call rnset(iseed)
                  
      do ir = 1, iiniter
          write(*,*) ir
          call gibbs_gj(iseed) 
          call gibbs_wi(iseed)          
          call gibbs_taui(iseed)          
      enddo
                         
      icount = 0
      do ir = 1, initer
          call gibbs(iseed)
      enddo    	  
                  
      do j = 1, nj
          do l = 1, nl
              cngj(j,l) = 0.d0
          enddo
      enddo                        
            
      do jj = 1, nc
	    mean_betaa(jj) = 0.d0          
      enddo
      do j = 1, nj
          do jj = 1, nb
	        mean_delta(j,jj) = 0.d0          
          enddo
      enddo
             
      bardic = 0.0d0      
      
      icount = 0
      do ir = 1, niter
              
          call gibbs(iseed)

          do l = 1, nl
              ngj(l) = 0
              do j = 1, nj
                  if (nint(gj(j)) .eq. l) then
                      ngj(l) = ngj(l) + 1
                  endif
              enddo
          enddo
                                                
          if (ir - ir/5*5 .eq. 1) then
        
              icount = icount + 1
              
              write(11,1) icount, ngj
              write(12,2) icount, betaa
              write(13,2) icount, alphaa, azeta, agamma
              write(14,2) icount, weight
              write(15,2) icount, theta
              write(16,2) icount, gj
                        
              do j = 1, nj
                  do l = 1, nl
                      if (l .eq. nint(gj(j))) then
                          cngj(j,l) = cngj(j,l) + 1.d0/dfloat(niter1)
                      endif
                  enddo
              enddo
                                      
              do jj = 1, nc
	            mean_betaa(jj) = mean_betaa(jj) 
     +                           + betaa(jj)/dfloat(niter1)
              enddo
              
              do j = 1, nj
                  
                  l = nint(gj(j))
                  
                  do jj = 1, nb
                      
                      delta(j,jj) = alphaa(jj,l)
                      
                      mean_delta(j,jj) = mean_delta(j,jj) 
     +                                 + delta(j,jj)/dfloat(niter1)
                      
                  enddo
                  
              enddo
            
              bardic = bardic 
     +               - 2.d0*gen_dic(betaa,delta)/dfloat(niter1)
              
          endif
                                   
      enddo    	  
      
      call rnget(iseed)
                                    
      do j = 1, nj
                                
          write(17,2) j, (cngj(j,l),l = 1, nl)
          
      enddo              
      
      dicbar = -2.d0*gen_dic(mean_betaa,mean_delta)
      pd = bardic - dicbar
      dic = dicbar + 2.0d0*pd 

      write(18,3) dicbar, pd, dic
                         
    1 format(1000i5)
    2 format(i5,1000f20.10)
    3 format(1000f20.10)
            
      end program
               
      
      real*8 function gen_dic(betaa,delta)
      implicit real*8 (a-h, o-z)
	integer, parameter :: ni = 11394, nj = 50
	integer, parameter :: nc = 5, nb = 3, nl = 2
		
	real*8 yi(ni), zi(ni,nc), xij(ni,nj,nb)    
      real*8 vnu, betaa(nc), delta(nj,nb)      
      
      real*8 zb, xa, fxij, cdf, dic
      
      common /vyi/yi
      common /vzi/zi
      common /vxij/xij
      
      common /vvnu/vnu
      
      external dtdf
                        
      dic = 0.d0
      do i = 1, ni
            
          zb = 0.d0
          do jj = 1, nc
              zb = zb + zi(i,jj)*betaa(jj)
          enddo
                    
          fxij = 0.d0               
          do j = 1, nj  
              if (xij(i,j,1) .ne. 0.d0) then
                                                  
                  xa = 0.d0
                  do jj = 1, nb
                      xa = xa + xij(i,j,jj)*delta(j,jj)
                  enddo
                      
                  fxij = fxij + xa
                  
              endif
          enddo
          
          cdf = dtdf(zb + fxij, vnu)
                  
          if (yi(i) .eq. 1.d0) then
                                            
              dic = dic + dlog(cdf)
                      
          else

              dic = dic + dlog(1.d0 - cdf)
                      
          endif
              
      enddo
                        
	gen_dic = dic

	end function          
      
 	include 'PesticideF_Gibbs.f'
	include 'tnorm.f'
	include 'optim1.f'
      include 'rGiGDist.f'                        
            