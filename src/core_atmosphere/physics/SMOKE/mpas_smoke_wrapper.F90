!>\file mpas_smoke_wrapper.F90
!! This file is MPAS Smoke wrapper
!! Haiqin.Li@noaa.gov 09/2024

module mpas_smoke_wrapper

   use mpas_kind_types
   use mpas_smoke_config
   use module_plumerise,      only : ebu_driver
   use module_add_emiss_burn, only : add_emis_burn
   use dep_dry_simple_mod,    only : dry_dep_driver_simple
   use dep_dry_emerson_mod,   only : dry_dep_driver_emerson
   use module_wetdep_ls,      only : wetdep_ls

   implicit none

   private

   public :: mpas_smoke_driver

contains

    subroutine mpas_smoke_driver(  smoke                 ,                              &
           ktau                  , dt                    , dxcell               ,       &
           xland                 , u10                   , v10                  ,       &
           ust                   , xlat                  , xlong                ,       &
           tskin                 , pblh                  , t2m                  ,       &
           p8w                   , dz8w                  , z_at_w               ,       &
           p_phy                 , t_phy                 , u_phy                ,       &
           v_phy                 , qv                    , vvel                 ,       &
           pi_phy                , rho_phy               , kpbl                 ,       &
           nsoil                 , smois                 , tslb                 ,       &
           ivgtyp                , isltyp                , nlcat                ,       &
           swdown                , z0                    , snowh                ,       &
           julian                , rmol                  , raincv               ,       &
           rainncv               , dpt2m                 , znt                  ,       &
           mavail                , g                     ,                              &
           cp                    , rd                    , gmt                  ,       &
           ids       , ide       , jds       , jde       , kds       , kde      ,       &
           ims       , ime       , jms       , jme       , kms       , kme      ,       &
           its       , ite       , jts       , jte       , kts       , kte              &
                            )
        
    implicit none

    !intent arguments:
    integer,intent(in):: ids,ide,jds,jde,kds,kde,        &
                         ims,ime,jms,jme,kms,kme,        &
                         its,ite,jts,jte,kts,kte
    integer,intent(in):: nsoil, nlcat, ktau
    integer,intent(in), dimension(ims:ime, jms:jme) :: isltyp, ivgtyp, kpbl

    real(RKIND),intent(in) :: dt, julian, g, cp, rd, gmt

    real(RKIND),intent(in), dimension(ims:ime, jms:jme) :: dxcell, xland, u10, v10,    &
                                        xlat, xlong, tskin, pblh, t2m, dpt2m
    real(RKIND),intent(in), dimension(ims:ime, jms:jme) :: swdown, z0, snowh, znt
    real(RKIND),intent(in), dimension(ims:ime, jms:jme) :: raincv, rainncv, mavail
                                        
    real(RKIND),intent(in), dimension(ims:ime, nsoil, jms:jme) :: smois, tslb
    real(RKIND),intent(in), dimension(ims:ime, kms:kme, jms:jme) :: p8w, dz8w, z_at_w,   &
                                        p_phy, t_phy, u_phy, v_phy, pi_phy, rho_phy, vvel
    real(RKIND),intent(inout), dimension(ims:ime, jms:jme) :: rmol, ust
    real(RKIND),intent(inout), dimension(ims:ime, kms:kme, jms:jme) :: qv
    real(RKIND),intent(inout), dimension(ims:ime, kms:kme, jms:jme) :: smoke


!>-- Local Variables
    real(RKIND), dimension(ims:ime, kms:kme, jms:jme, 1:num_chem )  :: chem
    real(RKIND), dimension(ims:ime, kms:kme, jms:jme) :: rri,     &
                     wind_phy,theta_phy,zmid,t8w
!>-- indexes, time
    integer :: julday
!>- dust & chemistry variables
    real(RKIND), dimension(ims:ime, jms:jme) :: ssm, rdrag, uthr  ! fengsha dust
    real(RKIND), dimension(ims:ime, nlcat, jms:jme) :: vegfrac
!>- plume variables
    ! -- buffers
    real(RKIND), dimension(ims:ime, kms:kme, jms:jme) :: ebu
    real(RKIND), dimension(ims:ime, kms:kme, jms:jme) :: ebu_smoke
    real(RKIND), dimension(ims:ime, jms:jme) :: ebu_in, frp_in, frp_inst
    real(RKIND), dimension(ims:ime, jms:jme) :: coef_bb_dc, flam_frac, frp,              &
                                                fire_hist, peak_hr, wetdpr_smoke,        &
                                                fire_end_hr, hwp_day_avg,                &
                                                uspdavg2d, hpbl2d, hwp_local
    real(RKIND), dimension(ims:ime, jms:jme) :: lu_nofire, lu_qfire, lu_sfire
    integer,     dimension(ims:ime, jms:jme) :: min_fplume2, max_fplume2, fire_type,     &
                                                kpbl_thetav
    logical :: call_plume
    real(RKIND), parameter :: conv_frpi   = 1.e-06_RKIND  ! FRP conversion factor, MW to W

!>- optical variables
    real(RKIND), dimension(ims:ime, jms:jme, ndvel) :: ddvel, settling_flux, drydep_flux_local
    real(RKIND), dimension(ims:ime, kms:kme, jms:jme, ndvel) :: vgrav

    real(RKIND) :: theta
    real(RKIND), dimension(ims:ime, kms:kme, jms:jme) :: thetav
!> -- other
    real(RKIND) :: curr_secs
    integer :: nbegin, nv
    integer :: i, j, k, kp, n

    character(100) :: errmsg
    integer        :: errflg

! FRP Thresholds
    real(RKIND), parameter :: frp_min        = 1.e+7     ! Minimum FRP (Watts) to distribute smoke in PBL, 10MW
    real(RKIND), parameter :: frp_max        = 2.e+10    ! Maximum FRP over 3km Pixel, 20,000 MW 
    real(RKIND), parameter :: zpbl_threshold = 2.e+3     ! Minimum PBL depth to have plume rise 
    real(RKIND), parameter :: uspd_threshold = 5.        ! Wind speed averaged across PBL depth to control smoke release levels 
    real(RKIND), parameter :: frp_wthreshold = 1.e+9     ! Minimum FRP (Watts) to have plume rise in windy conditions
    real(RKIND), parameter :: ebb_min        = 1.e-3     ! Minimum smoke emissions (ug/m2/s)

    !read(pram.xt,...)


    errmsg = ''
    errflg = 0

    if (.not. do_mpas_smoke) return

    min_fplume2 = 0
    max_fplume2 = 0
    uspdavg2d   = 0.
    hpbl2d      = 0.
    peak_hr     = 0.
    fire_type   = 0
    flam_frac   = 0.


    curr_secs = ktau * dt
    julday = int(julian)

    ! plumerise frequency in minutes set up by the namelist input
    call_plume       = (do_plumerise .and. (plumerisefire_frq > 0))
    if (call_plume) call_plume = (mod(int(curr_secs), max(1, 60*plumerisefire_frq)) == 0) .or. (ktau == 2)
    
    do j=jts,jte
    do i=its,ite
     if (xland(i,j).eq.1.)then
      ebu_in(i,j) = 100. !smoke_RRFS(i,hour_int+1,1) ! smoke
      frp_in(i,j) = 100. !smoke_RRFS(i,hour_int+1,2)*conv_frp ! frp
     endif
    enddo
    enddo

    !>- get ready for chemistry run
    call mpas_smoke_prep(                                                   &
        ktau, nlcat,cp,ebb_dcycle,ebb_min,                                  &
        xland,xlat,xlong,ivgtyp,isltyp,vegfrac,                             &
        snowh,u10,v10,t2m,dpt2m,mavail,hwp_local,hwp_day_avg,               &
        ebu_in,t_phy,u_phy,v_phy,p_phy,pi_phy,z_at_w,qv,                    &
        wind_phy,theta_phy,zmid,kpbl_thetav,                                &
        peak_hr,coef_bb_dc,fire_hist,                                       &
        lu_nofire, lu_qfire, lu_sfire, fire_type,                           &
        smoke,num_chem,chem,                                                &
        ids,ide, jds,jde, kds,kde,                                          &
        ims,ime, jms,jme, kms,kme,                                          &
        its,ite, jts,jte, kts,kte)

    if (ktau==1) then
      ebu = 0.
      do j=jts,jte
      do i=its,ite
      do k=kts+1,kte
         ebu(i,k,j)= 0.
      enddo
      enddo
      enddo
    else
      do j=jts,jte
      do k=kts,kte
      do i=its,ite
      ! ebu is divided by coef_bb_dc since it is applied in the output
        ebu(i,k,j)=ebu_smoke(i,k,j) / MAX(1.E-4,coef_bb_dc(i,j))
      enddo
      enddo
      enddo
    endif
   
    

    ! compute wild-fire plumes
    if (call_plume) then
      ! Apply the diurnal cycle coefficient to frp_inst ()
      do j=jts,jte
      do i=its,ite
        if ( fire_type(i,j) .eq. 4 ) then ! only apply scaling factor to wildfires
           frp_inst(i,j) = min(sc_factor*frp_in(i,j)*coef_bb_dc(i,j),frp_max)
        else
           frp_inst(i,j) = min(frp_in(i,j)*coef_bb_dc(i,j),frp_max)
        endif
      enddo
      enddo

      call ebu_driver (                                               &
                 flam_frac,ebu_in,ebu,                                &
                 theta_phy,qv,                                        &
                 rho_phy,vvel,u_phy,v_phy,pi_phy,wind_phy,            &
                 z_at_w,zmid,g,cp,rd,                                 &
                 frp_inst, min_fplume2, max_fplume2,                  &
                 plume_wind_eff,                                      &
                 kpbl_thetav,kpbl,curr_secs,                          &
                 xlat, xlong, uspdavg2d, hpbl2d, plume_alpha,         &
                 frp_min, frp_wthreshold,                             &
                 zpbl_threshold, uspd_threshold,                      &
                 ids,ide, jds,jde, kds,kde,                           &
                 ims,ime, jms,jme, kms,kme,                           &
                 its,ite, jts,jte, kts,kte, errmsg, errflg            )
      if(errflg/=0) return
    end if

    ! -- add biomass burning emissions at every timestep
    if (addsmoke_flag == 1) then
     call add_emis_burn(dt,dz8w,rho_phy,pi,ebb_min,                   &
                        chem,julday,gmt,xlat,xlong,                   &
                        fire_end_hr, peak_hr,curr_secs,               &
                        coef_bb_dc,fire_hist,hwp_local,hwp_day_avg,   &
                        swdown,ebb_dcycle,ebu_in,ebu,fire_type,       &
                        qv, add_fire_moist_flux,                      &
                        sc_factor,                                    &    
                        ids,ide, jds,jde, kds,kde,                    &
                        ims,ime, jms,jme, kms,kme,                    &
                        its,ite, jts,jte, kts,kte                     )
    endif


    !>-- compute dry deposition, based on Emerson et al., (2020)
    if (drydep_opt == 1) then
     call dry_dep_driver_emerson(rmol,ust,znt,ndvel,ddvel,            &
        vgrav,chem,dz8w,snowh,t_phy,p_phy,rho_phy,ivgtyp,g,dt,        &
        drydep_flux_local,settling_flux,dbg_opt,                      &
        ids,ide, jds,jde, kds,kde,                                    &
        ims,ime, jms,jme, kms,kme,                                    &
        its,ite, jts,jte, kts,kte, curr_secs, xlat, xlong      )
    !>-- compute dry deposition based on simple parameterization (HRRR-Smoke)
    elseif (drydep_opt == 2) then
     call dry_dep_driver_simple(rmol,ust,ndvel,ddvel,                 &
        ids,ide, jds,jde, kds,kde,                                    &
        ims,ime, jms,jme, kms,kme,                                    &
        its,ite, jts,jte, kts,kte                                     )
    else
        ddvel=0.
    endif

    !>- large-scale wet deposition
    if (wetdep_ls_opt == 1) then
       call  wetdep_ls(dt,chem,rainncv,qv,                            &
                     rho_phy,num_chem,ndvel,dz8w,vvel,                &
                     wetdpr_smoke,                                    &
                     ids,ide, jds,jde, kds,kde,                       &
                     ims,ime, jms,jme, kms,kme,                       &
                     its,ite, jts,jte, kts,kte                        )
    endif

    !>-- output of MPAS-Smoke
    do j=jts,jte
    do k=kts,kte
    do i=its,ite
       ebu_smoke(i,k,j)=ebu(i,k,j) * coef_bb_dc(i,j)
    enddo
    enddo
    enddo

!---- put smoke stuff back into tracer array
    do j=jts,jte
    do k=kts,kte
    do i=its,ite
      smoke(i,k,j)  = min(5000.,max(epsilc,chem(i,k,j,p_smoke))) 
    enddo
    enddo
    enddo

 end subroutine mpas_smoke_driver

 subroutine mpas_smoke_prep(                                                &
        ktau, nlcat,cp,ebb_dcycle,ebb_min,                                  &
        xland,xlat,xlong,ivgtyp,isltyp,vegfrac,                             &
        snowh,u10,v10,t2m,dpt2m,wetness,hwp_local,hwp_day_avg,              &
        ebu_in,t_phy,u_phy,v_phy,p_phy,pi_phy,z_at_w,qv,                    &
        wind_phy,theta_phy,zmid,kpbl_thetav,                                &
        peak_hr,coef_bb_dc,fire_hist,                                       &
        lu_nofire, lu_qfire, lu_sfire, fire_type,                           &
        smoke,num_chem,chem,                                                &
        ids,ide, jds,jde, kds,kde,                                          &
        ims,ime, jms,jme, kms,kme,                                          &
        its,ite, jts,jte, kts,kte)

    !intent arguments:
    integer,intent(in):: ids,ide,jds,jde,kds,kde,                           &
                         ims,ime,jms,jme,kms,kme,                           &
                         its,ite,jts,jte,kts,kte

    integer,intent(in):: ktau, nlcat,ebb_dcycle,num_chem
    integer,intent(in), dimension(ims:ime, jms:jme) :: isltyp, ivgtyp

    real(RKIND),intent(in) :: cp, ebb_min
    real(RKIND), dimension(ims:ime, nlcat, jms:jme) :: vegfrac

    real(RKIND),intent(in), dimension(ims:ime, jms:jme) :: xland, xlat, xlong, ebu_in,           &
                                        snowh, u10, v10, t2m, dpt2m, wetness
    real(RKIND),intent(in), dimension(ims:ime, kms:kme, jms:jme) :: qv, z_at_w,                  &
                                        p_phy, t_phy, u_phy, v_phy, pi_phy
    real(RKIND),intent(out), dimension(ims:ime, kms:kme, jms:jme) :: zmid, wind_phy,theta_phy

    integer    ,intent(out), dimension(ims:ime, jms:jme) :: fire_type, kpbl_thetav
    real(RKIND),intent(out), dimension(ims:ime, jms:jme) :: lu_nofire, lu_qfire, lu_sfire
    real(RKIND),intent(out), dimension(ims:ime, jms:jme) :: coef_bb_dc, fire_hist, peak_hr,      &
                                                            hwp_local,hwp_day_avg
    real(RKIND),intent(inout), dimension(ims:ime, kms:kme, jms:jme) :: smoke
    real(RKIND),intent(inout), dimension(ims:ime, kms:kme, jms:jme, 1:num_chem )  :: chem

    !local variables
    real(RKIND), parameter :: delta_theta4gust = 0.5
    real(RKIND) :: theta, wdgust, snoweq
!    real(RKIND), dimension(ims:ime, jms:jme) :: 
    real(RKIND), dimension(ims:ime, kms:kme, jms:jme) :: thetav

    integer :: i, j, k, k1

    !if (ktau==1) then
      do j=jts,jte
      do i=its,ite
        fire_hist   (i,j) = 1.
        coef_bb_dc  (i,j) = 1.
        if (xlong(i,j)<230.) then
            peak_hr(i,j)= 0.0* 3600.     ! peak at 24 UTC, fires in Alaska
        elseif(xlong(i,j)<245.) then
            peak_hr(i,j)= 23.0* 3600.
        elseif (xlong(i,j)<260.) then
            peak_hr(i,j)= 22.0* 3600.    ! peak at 22 UTC, fires in the western US
        elseif (xlong(i,j)<275.) then
            peak_hr(i,j)= 21.0* 3600.
        elseif (xlong(i,j)<290.) then    ! peak at 20 UTC, fires in the eastern US
            peak_hr(i,j)= 20.0* 3600.
        else
            peak_hr(i,j)= 19.0* 3600.
        endif
      enddo
      enddo
    !endif

    do j=jts,jte
    do k=kts,kte
    do i=its,ite
      zmid(i,k,j)=z_at_w(i,k,j)
    enddo
    enddo
    enddo

    do j=jts,jte
    do k=kts,kte
    do i=its,ite
      theta_phy(i,k,j) = t_phy(i,k,j)/pi_phy(i,k,j)*cp
      wind_phy(i,k,j) = sqrt(u_phy(i,k,j)**2 + v_phy(i,k,j)**2)
    enddo
    enddo
    enddo

    !---- Calculate PBLH and K-PBL based on virtual potential temperature profile
    !---- First calculate THETAV
    do j = jts,jte
    do k = kts,kte
    do i = its,ite
       theta = t_phy(i,k,j) * (1.E5/p_phy(i,k,j))**0.286
       thetav(i,k,j) = theta * (1. + 0.61 * qv(i,k,j))
    enddo
    enddo
    enddo
    !---- Now use the UPP code to deterimine the height and level
    do i = its, ite
    do j = jts, jte
       if ( thetav(i,kts+1,j) .lt. ( thetav(i,kts,j) + delta_theta4gust) ) then
          do k = kts+1, kte
             k1 = k
!--- give theta-v at the sfc a 0.5K boost in the PBLH definition
             if ( thetav(i,kts+k-1,j) .gt. ( thetav(i,kts,j) + delta_theta4gust) ) then
                exit
             endif
          enddo
          kpbl_thetav(i,j) = k1
       else
          kpbl_thetav(i,j) = kts + 1
       endif
   enddo
   enddo

    !RAR: change this to the fractional LU type; fire_type: 0- no fires, 1- Ag
! or urban fires, 2- prescribed fires in wooded area, 3- wildfires
    if (ebb_dcycle==2) then
      do j=jts,jte
      do i=its,ite
        if (ebu_in(i,j)<ebb_min) then
           fire_type(i,j) = 0
           lu_nofire(i,j) = 1.0
        else
          ! Permanent wetlands, snow/ice, water, barren tundra:
          lu_nofire(i,j)= vegfrac(i,11,j) + vegfrac(i,15,j) + vegfrac(i,17,j) + vegfrac(i,20,j)
          ! cropland, urban, cropland/natural mosaic, barren and sparsely
          ! vegetated and non-vegetation areas:
          lu_qfire(i,j) = lu_nofire(i,j) + vegfrac(i,12,j) + vegfrac(i,13,j) + vegfrac(i,14,j) + vegfrac(i,16,j)
          ! Savannas and grassland fires, these fires last longer than the Ag fires:
          lu_sfire(i,j) = lu_qfire(i,j) + vegfrac(i,8,j) + vegfrac(i,9,j) + vegfrac(i,10,j)
          if (lu_nofire(i,j)>0.95) then ! no fires
            fire_type(i,j) = 0
          else if (lu_qfire(i,j)>0.9) then   ! Ag. and urban fires
            fire_type(i,j) = 1
          else if (xlong(i,j)>260. .AND. xlat(i,j)>25. .AND. xlat(i,j)<41.) then
            fire_type(i,j) = 2    ! slash burn and wildfires in the east, eastern temperate forest ecosystem
          else if (lu_sfire(i,j)>0.8) then
            fire_type(i,j) = 3    ! savanna and grassland fires
          else
            fire_type(i,j) = 4    ! potential wildfires
          end if
        end if
      end do
      end do
    endif

    !>-- HWP: Pre-release of RRFSv1 method - using wind gust calculated via UPP Method
    do i=its, ite
    do j=jts, jte
      wdgust  =max(sqrt(u10(i,j)**2+v10(i,j)**2),3.)
      snoweq  =max((25.-snowh(i,j))/25.,0.)
      hwp_local(i,j)=0.177*wdgust**0.97*max(t2m(i,j)-dpt2m(i,j),15.)**1.03*((1.-wetness(i,j))**0.4)*snoweq 
      hwp_day_avg(i,j)=hwp_local(i,j)
    enddo
    enddo


    do j=jms,jte
    do k=kms,kte
    do i=ims,ime
      chem(i,k,j,p_smoke)=max(epsilc,smoke(i,k,j))
    enddo
    enddo
    enddo

  end subroutine mpas_smoke_prep

!> @}
  end module mpas_smoke_wrapper
