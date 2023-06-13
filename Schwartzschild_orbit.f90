!!!******************************************************************************
!!!*						                	                                *
!!!*                         Schwartzschild Orbit      		                    *
!!!*							        	                                    *
!!!*                        Author: Dr. Andrea Boghi 		        	        *
!!!*							        	                                    *
!!!*	(dv/dOmega)**2 -delta*v**3 +v**2 -2*v*alpha = e**2 -1 + alpha*2/delta	*
!!!*							        	                                    *
!!!*                    d2v/dOmega2 + v -1.5*delta*v**2 = alpha	                *
!!!*								                                            *
!!!******************************************************************************
    
Program Schwartzschild_Orbit
implicit none
real*8 Omega,v(1:4),v0(1:4),v_temp(1:4)
real*8 dvdOmega_1(1:4),dvdOmega_2(1:4),dvdOmega_3(1:4),dvdOmega_4(1:4)
integer i,j,k,n
real*8  dOmega,delta,e,PI,alpha,r
logical is_particle
character(len=80) :: file_name

open(1,file='imp.dat')
read(1,*) is_particle
read(1,*) k
read(1,*) dOmega
read(1,*) delta
read(1,*) e
close(1)

print*, ' is_particle =', is_particle, ' number of cycles =', k
print*, ' delta =', delta, ' e =', e, ' dOmega [deg] =', dOmega

n = INT(k * 360 / dOmega)

PI = 2.d0 * ASIN(1.d0)

dOmega = dOmega * 2.d0 * PI / 360.d0

print*, ' PI =', PI, ' n =', n, ' dOmega =', dOmega

Omega = 0.d0
v(1) = 1.d0
if(is_particle) then
    alpha = 1.d0
    write(file_name,121)'Particle_e_',e,'_delta_',delta,'.csv'
    121 format(a11,e9.3,a7,e9.3,a4)
else 
    alpha = 0.d0
    write(file_name,122)'Photon_e_',e,'_delta_',delta,'.csv'
    122 format(a9,e9.3,a7,e9.3,a4)
endif

v(2) = ( e**2 + delta * v(1)**3.d0 - delta - v(1)**2.d0 + 1.d0 + 2.d0 * alpha * ( v(1) - 1.d0) )**0.5d0 
v(3) = 0.d0
v(4) = 0.d0

open(11,file=file_name,form='formatted')
write(11,*) '"Omega","r","x","y","tau","t"'
write(11,101) Omega,',',1.d0/v(1),',',(1.d0/v(1))*cos(Omega),',',(1.d0/v(1))*sin(Omega),',',v(3),',',v(4)

print*, ' v(2) = ', v(2)

j=1
do while(j<=n .and. ( abs(v(1)) < 1.d+30 .and. v(1) > 1.d-03 ) )
    
    do i=1,4
       v0(i) = v(i)
    enddo
    Omega = Omega + dOmega

    !!! First Step
    call Schwartzschild_Model(dvdOmega_1,v0,delta,alpha,e)
    
    do i=1,4
       v_temp(i) = v0(i) + ( dOmega / 2.d0 ) * dvdOmega_1(i)
    enddo

    !!! Second Step
    call Schwartzschild_Model(dvdOmega_2,v_temp,delta,alpha,e)
      
    do i=1,4
       v_temp(i) = v0(i) + ( dOmega / 2.d0 ) * dvdOmega_2(i)
    enddo

    !!! Third Step
    call Schwartzschild_Model(dvdOmega_3,v_temp,delta,alpha,e)

    do i=1,4
       v_temp(i) = v0(i) + dOmega * dvdOmega_3(i)
    enddo
    
    !!! Fourth Step
    call Schwartzschild_Model(dvdOmega_4,v_temp,delta,alpha,e)

    do i=1,4
       v(i) = v0(i) + ( dOmega / 6.d0 ) * ( dvdOmega_1(i) + 2.d0 * dvdOmega_2(i) + 2.d0 * dvdOmega_3(i) + dvdOmega_4(i) )
    enddo
   
    if(v(1)>0.d0) then
      write(11,101) Omega,',',1.d0/v(1),',',(1.d0/v(1))*cos(Omega),',',(1.d0/v(1))*sin(Omega),',',v(3),',',v(4)
    endif
    print*, ' j =', j,' v(1) =', v(1),' v(2) =', v(2),' v(3) =', v(3),' v(4) =', v(4)
    j=j+1

enddo

close(11)

101 format(e18.10,A,e18.10,A,e18.10,A,e18.10,A,e18.10,A,e18.10)

end

!!!******************************************************************************
!!!*						                	                                *
!!!*                         Schwartzschild Model     		                    *
!!!*							        	                                    *
!!!*	(dv/dOmega)**2 -delta*v**3 +v**2 -2*v*alpha = e**2 -1 + alpha*2/delta	*
!!!*							        	                                    *
!!!*                    d2v/dOmega2 + v -1.5*delta*v**2 = alpha	                *
!!!*								                                            *
!!!******************************************************************************
      
subroutine  Schwartzschild_Model(dvdOmega,v,delta,alpha,e)
implicit none
real*8, intent(in) :: v(1:4),delta,alpha,e
real*8, intent(out) :: dvdOmega(1:4)

dvdOmega(1) = v(2)
dvdOmega(2) = -v(1) + 1.5d0 * delta * v(1)**2.d0 +alpha 
dvdOmega(3) = 1.d0/v(1)**2.d0
dvdOmega(4) = ( ( alpha * (1.d0 - delta) + (delta/2.d0) * (e**2.d0 + 1.d0 - delta) )**5.d-1 ) &
/( ( 1.d0 - delta * v(1) ) * v(1)**2.d0 )

end
