!>\file module_wetdep_ls.F90
!! This file contains aerosol wet deposition module.

module module_wetdep_ls
  use mpas_kind_types
  use mpas_smoke_config, only : p_smoke, alpha => wetdep_ls_alpha

contains
subroutine wetdep_ls(dt,var,rain,moist,                                         &
                     rho,nchem,ndvel,dz8w,vvel,                                 &
                     wetdpr_smoke,                                              &
                     ids,ide, jds,jde, kds,kde,                                 &
                     ims,ime, jms,jme, kms,kme,                                 &
                     its,ite, jts,jte, kts,kte )
   implicit none

   integer,      intent(in) :: nchem, ndvel,                                    &
                               ids,ide, jds,jde, kds,kde,                       &
                               ims,ime, jms,jme, kms,kme,                       &
                               its,ite, jts,jte, kts,kte
   real(RKIND), intent(in) :: dt
   real(RKIND), dimension( ims:ime, kms:kme, jms:jme),intent(in) :: moist
   real(RKIND), dimension( ims:ime, kms:kme, jms:jme),intent(in) :: rho,dz8w,vvel        
   real(RKIND), dimension( ims:ime, kms:kme, jms:jme,1:nchem),intent(inout) :: var        
   real(RKIND), dimension( ims:ime, jms:jme ), intent(out) :: wetdpr_smoke
   real(RKIND), dimension( ims:ime, jms:jme),intent(in) :: rain
   real(RKIND), dimension( its:ite, jts:jte) :: var_sum,var_rmv
   real(RKIND), dimension( its:ite, kts:kte, jts:jte) :: var_rmvl
   real(RKIND), dimension( its:ite, jts:jte) :: frc,var_sum_clw,rain_clw     
   real(RKIND) :: dvar,factor,clsum
   integer :: nv,i,j,k,km,kb,kbeg
  !real(RKIND), parameter :: alpha = .5 ! scavenging factor
   integer, save :: print_alpha = 0

    wetdpr_smoke   =0.

    !if ( print_alpha == 0 ) then
    !   write(*,*) 'wetdep_ls, alpha = ',alpha
    !   print_alpha = print_alpha + 1
    !endif


    do nv=1,nchem
      do i=its,ite
       do j=jts,jte
        var_sum_clw(i,j)=0.
        var_sum(i,j)=0.
        var_rmvl(i,:,j)=0.
        frc(i,j)=0.
        rain_clw(i,j)=0.
        if(rain(i,j).gt.1.e-10)then
! convert rain back to rate
!
           rain_clw(i,j)=rain(i,j)/dt
! total cloud water
!
           do k=1,kte-1
              dvar=max(0.,moist(i,k,j)*rho(i,k,j)*vvel(i,k,j)*dz8w(i,k,j))
              var_sum_clw(i,j)=var_sum_clw(i,j)+dvar
              var_sum(i,j)=var_sum(i,j)+var(i,k,j,nv)*rho(i,k,j)
           enddo
           if(var_sum(i,j).gt.1.e-10 .and. var_sum_clw(i,j).gt.1.e-10 ) then
!          assuming that frc is onstant, it is my conversion factor 
!          (just like in convec. parameterization)
              frc(i,j)=rain_clw(i,j)/var_sum_clw(i,j)
              frc(i,j)=max(1.e-6,min(frc(i,j),.005))
           endif
        endif
       enddo
       enddo
!
! get rid of it
!
       do i=its,ite
       do j=jts,jte
       if(rain(i,j).gt.1.e-10 .and. var_sum(i,j).gt.1.e-10 .and. var_sum_clw(i,j).gt.1.e-10)then
         do k=kts,kte-2
          if(var(i,k,j,nv).gt.1.e-16 .and. moist(i,k,j).gt.0.)then
            factor = max(0.,frc(i,j)*rho(i,k,j)*dz8w(i,k,j)*vvel(i,k,j))
            dvar=alpha*factor/(1+factor)*var(i,k,j,nv)
! Accumulate diags
            if (nv .eq. p_smoke ) then
               wetdpr_smoke(i,j) = wetdpr_smoke(i,j) + dvar * rho(i,k,j) / dt
            !elseif (nv .eq. p_dust_1 ) then
            !   wetdpr_dust(i,j) = wetdpr_dust(i,j) + dvar * rho(i,k,j) / dt
            !elseif (nv .eq. p_coarse_pm ) then
            !   wetdpr_coarsepm(i,j) = wetdpr_coarsepm(i,j) + dvar * rho(i,k,j) / dt
            endif
            var(i,k,j,nv)=max(1.e-16,var(i,k,j,nv)-dvar)
          endif
         enddo
       endif
       enddo
       enddo
      enddo ! nv
end subroutine wetdep_ls 
end module module_wetdep_ls
