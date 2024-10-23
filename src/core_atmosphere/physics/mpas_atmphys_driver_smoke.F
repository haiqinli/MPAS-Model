! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!=================================================================================================================
 module mpas_atmphys_driver_smoke
 use mpas_kind_types
 use mpas_pool_routines
 use mpas_timer, only : mpas_timer_start, mpas_timer_stop
 use mpas_atmphys_manager, only: gmt,curr_julday,julday,year

 use mpas_atmphys_constants
 use mpas_atmphys_utilities
 use mpas_atmphys_vars

 use mpas_smoke_wrapper

 implicit none
 private
 public:: allocate_smoke,     &
          deallocate_smoke,   &
          init_smoke,         &
          driver_smoke

!MPAS driver for smoke.
!Haiqin Li (haiqin.li@noaa.gov).
!2024-09-25.

 contains

!=================================================================================================================
 subroutine allocate_smoke(configs)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: configs

 integer:: i,k,j

 if(.not.allocated(dx_p)        ) allocate(dx_p(ims:ime,jms:jme)                 )
 if(.not.allocated(xland_p)     ) allocate(xland_p(ims:ime,jms:jme)              )
 if(.not.allocated(u10_p)       ) allocate(u10_p(ims:ime,jms:jme)                )
 if(.not.allocated(v10_p)       ) allocate(v10_p(ims:ime,jms:jme)                )
 if(.not.allocated(ust_p)       ) allocate(ust_p(ims:ime,jms:jme)                )
 if(.not.allocated(xlat_p)      ) allocate(xlat_p(ims:ime,jms:jme)               )
 if(.not.allocated(xlon_p)      ) allocate(xlon_p(ims:ime,jms:jme)               )
 if(.not.allocated(vegfra_p)    ) allocate(vegfra_p(ims:ime,jms:jme)             )

 if(.not.allocated(tsk_p)       ) allocate(tsk_p(ims:ime,jms:jme)                )
 if(.not.allocated(hpbl_p)      ) allocate(hpbl_p(ims:ime,jms:jme)               )
 if(.not.allocated(t2m_p)       ) allocate(t2m_p(ims:ime,jms:jme)                )
 if(.not.allocated(qv_p)        ) allocate(qv_p(ims:ime,kms:kme,jms:jme)         )
 if(.not.allocated(w_p)         ) allocate(w_p(ims:ime,kms:kme,jms:jme)          )
 if(.not.allocated(isltyp_p)    ) allocate(isltyp_p(ims:ime,jms:jme)             )
 if(.not.allocated(swdown_p)    ) allocate(swdown_p(ims:ime,jms:jme)             )
 if(.not.allocated(z0_p)        ) allocate(z0_p(ims:ime,jms:jme)                 )
 if(.not.allocated(snowh_p)     ) allocate(snowh_p(ims:ime,jms:jme)              )
 if(.not.allocated(rmol_p)      ) allocate(rmol_p(ims:ime,jms:jme)               )
 if(.not.allocated(raincv_p)    ) allocate(raincv_p(ims:ime,jms:jme)             )
 if(.not.allocated(rainncv_p)   ) allocate(rainncv_p(ims:ime,jms:jme)            )
 if(.not.allocated(kpbl_p)      ) allocate(kpbl_p(ims:ime,jms:jme)               )
 if(.not.allocated(znt_p)       ) allocate(znt_p(ims:ime,jms:jme)                )
 if(.not.allocated(dpt2m_p)     ) allocate(dpt2m_p(ims:ime,jms:jme)              )
 if(.not.allocated(mavail_p)    )  allocate(mavail_p(ims:ime,jms:jme)            )

 if(.not.allocated(smois_p)     ) allocate(smois_p(ims:ime,1:num_soils,jms:jme)  )
 if(.not.allocated(tslb_p)      ) allocate(tslb_p(ims:ime,1:num_soils,jms:jme)   )
 if(.not.allocated(smoke_p)     ) allocate(smoke_p(ims:ime,kms:kme,jms:jme)      )

 end subroutine allocate_smoke

!=================================================================================================================
 subroutine deallocate_smoke(configs)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: configs

 if(allocated(dx_p)           ) deallocate(dx_p           )
 if(allocated(xland_p)        ) deallocate(xland_p        )
 if(allocated(u10_p)          ) deallocate(u10_p          )
 if(allocated(v10_p)          ) deallocate(v10_p          )
 if(allocated(ust_p)          ) deallocate(ust_p          )
 if(allocated(xlat_p)         ) deallocate(xlat_p         )
 if(allocated(xlon_p)         ) deallocate(xlon_p         )
 if(allocated(vegfra_p)       ) deallocate(vegfra_p       )

 if(allocated(tsk_p)          ) deallocate(tsk_p          )
 if(allocated(hpbl_p)         ) deallocate(hpbl_p         )
 if(allocated(t2m_p)          ) deallocate(t2m_p          )
 if(allocated(qv_p)           ) deallocate(qv_p           )
 if(allocated(w_p)            ) deallocate(w_p            )
 if(allocated(isltyp_p)       ) deallocate(isltyp_p       )
 if(allocated(swdown_p)       ) deallocate(swdown_p       )
 if(allocated(z0_p)           ) deallocate(z0_p           )
 if(allocated(snowh_p)        ) deallocate(snowh_p        )
 if(allocated(rmol_p)         ) deallocate(rmol_p         )
 if(allocated(raincv_p)       ) deallocate(raincv_p       )
 if(allocated(rainncv_p)      ) deallocate(rainncv_p      )
 if(allocated(kpbl_p)         ) deallocate(kpbl_p         )
 if(allocated(znt_p)          ) deallocate(znt_p          )
 if(allocated(dpt2m_p)        ) deallocate(dpt2m_p        )
 if(allocated(mavail_p)       ) deallocate(mavail_p       )

 if(allocated(smois_p)        ) deallocate(smois_p        )
 if(allocated(tslb_p)         ) deallocate(tslb_p         )
 if(allocated(smoke_p)        ) deallocate(smoke_p        )

!-----------------------------------------------------------------------------------------------------------------

 end subroutine deallocate_smoke

 !=================================================================================================================
 subroutine smoke_from_MPAS(dt_dyn,configs,mesh,sfc_input,diag_physics,tend_physics,its,ite)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: configs
 type(mpas_pool_type),intent(in):: mesh
 type(mpas_pool_type),intent(in):: diag_physics
 type(mpas_pool_type),intent(in):: sfc_input
 type(mpas_pool_type),intent(in):: tend_physics

 integer,intent(in):: its,ite
 real(kind=RKIND),intent(in):: dt_dyn

 !local pointers:
 real(kind=RKIND),dimension(:),pointer  :: meshDensity
 real(kind=RKIND),dimension(:),pointer  :: xland
 real(kind=RKIND),pointer               :: len_disp
 real(kind=RKIND),dimension(:),pointer  :: ust
 real(kind=RKIND),dimension(:),pointer  :: latCell,lonCell
 real(kind=RKIND),dimension(:),pointer  :: raincv,rainncv,mavail
 real(kind=RKIND),dimension(:),pointer  :: rmol,snowh,z0,t2m,dpt2m
 real(kind=RKIND),dimension(:),pointer  :: u10,v10,znt,hpbl,skintemp
 real(kind=RKIND),dimension(:,:),pointer:: smois,tslb
 real(kind=RKIND),dimension(:,:),pointer:: zgrid
 integer,dimension(:),pointer:: isltyp,ivgtyp,kpbl
 integer:: i,j,k,n

 call mpas_pool_get_config(configs,'config_len_disp',len_disp)
 call mpas_pool_get_array(mesh,'meshDensity',meshDensity)
 call mpas_pool_get_array(mesh,'latCell',latCell)
 call mpas_pool_get_array(mesh,'lonCell',lonCell)
 call mpas_pool_get_array(sfc_input,'xland',xland)
 call mpas_pool_get_array(diag_physics,'ust',ust)

 call mpas_pool_get_array(diag_physics,'u10',u10)
 call mpas_pool_get_array(diag_physics,'v10',v10)
 call mpas_pool_get_array(sfc_input,'skintemp',skintemp)
 call mpas_pool_get_array(diag_physics,'hpbl',hpbl)
 call mpas_pool_get_array(diag_physics,'kpbl',kpbl)
 call mpas_pool_get_array(diag_physics,'znt',znt)

 call mpas_pool_get_array(mesh,'zgrid',zgrid)
 call mpas_pool_get_array(sfc_input,'smois',smois)
 call mpas_pool_get_array(sfc_input,'tslb',tslb)
 call mpas_pool_get_array(sfc_input,'isltyp',isltyp)
 call mpas_pool_get_array(sfc_input,'ivgtyp',ivgtyp)
 call mpas_pool_get_array(diag_physics, 'z0',z0)
 call mpas_pool_get_array(sfc_input,'snowh',snowh)
 call mpas_pool_get_array(diag_physics,'rmol',rmol)
 call mpas_pool_get_array(diag_physics,'raincv',raincv)
 call mpas_pool_get_array(diag_physics,'rainncv',rainncv)
 call mpas_pool_get_array(diag_physics,'dewpoint_surface',dpt2m)
 call mpas_pool_get_array(diag_physics,'mavail',mavail)

 do j = jts,jte
 do n = 1,num_soils
 do i = its,ite
    smois_p(i,n,j)  = smois(n,i)
    tslb_p(i,n,j)   = tslb(n,i)
 enddo
 enddo
 enddo

 do j = jts,jte
 do i = its,ite
    dx_p(i,j)      = len_disp / meshDensity(i)**0.25
    xland_p(i,j)   = xland(i)
    ust_p(i,j)     = ust(i)
    xlat_p(i,j)    = latCell(i) / degrad
    xlon_p(i,j)    = lonCell(i) / degrad
    isltyp_p(i,j)  = isltyp(i)
    ivgtyp_p(i,j)  = ivgtyp(i)
    kpbl_p(i,j)    = kpbl(i)

    u10_p(i,j)     = u10(i)
    v10_p(i,j)     = v10(i)
    znt_p(i,j)     = znt(i)
    tsk_p(i,j)     = skintemp(i)
    hpbl_p(i,j)    = hpbl(i)
    z0_p(i,j)      = z0(i)
    snowh_p(i,j)   = snowh(i)
    rmol_p(i,j)    = rmol(i)
    raincv_p(i,j)  = raincv(i)
    rainncv_p(i,j) = rainncv(i)
    dpt2m_p(i,j)   = dpt2m(i)
    mavail_p(i,j)  = mavail(i)
    ! smoke_p
 enddo
 enddo

 end subroutine smoke_from_MPAS

!=================================================================================================================
 subroutine smoke_to_MPAS(configs,diag_physics,tend_physics,its,ite)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: configs

!inout arguments:
 type(mpas_pool_type),intent(inout):: diag_physics
 type(mpas_pool_type),intent(inout):: tend_physics

 integer,intent(in):: its,ite

!local variables:
 integer:: i,k,j

!local pointers:

!inout variables:


 end subroutine smoke_to_MPAS

!=================================================================================================================
 subroutine init_smoke(mesh,configs,diag_physics)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: mesh
 type(mpas_pool_type),intent(in):: configs

!inout arguments:
 type(mpas_pool_type),intent(inout):: diag_physics

 end subroutine init_smoke

!=================================================================================================================
 subroutine driver_smoke(itimestep,configs,mesh,sfc_input,diag_physics,tend_physics,its,ite)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: configs
 type(mpas_pool_type),intent(in):: mesh

 integer,intent(in):: its,ite
 integer,intent(in):: itimestep

!inout arguments:
 type(mpas_pool_type),intent(inout):: sfc_input
 type(mpas_pool_type),intent(inout):: diag_physics
 type(mpas_pool_type),intent(inout):: tend_physics

 integer:: ktau,initflag
 real(kind=RKIND):: curr_secs, dt

 !local pointers:
 logical,pointer:: config_do_restart
 real(kind=RKIND),dimension(:,:),pointer:: zgrid

!copy MPAS arrays to local arrays:
 call smoke_from_MPAS(dt_dyn,configs,mesh,sfc_input,diag_physics,tend_physics,its,ite)
 call mpas_pool_get_config(configs,'config_do_restart'       ,config_do_restart)
 call mpas_pool_get_array(mesh        ,'zgrid' ,zgrid)

 initflag = 1
 if(config_do_restart .or. itimestep > 1) initflag = 0



 call mpas_timer_start('mpas_smoke')
 call mpas_smoke_driver(           smoke    = smoke_p    ,                               &
           ktau    = itimestep   , dt       = dt_dyn     , dxcell   = dx_p       ,       &
           xland   = xland_p     , u10      = u10_p      , v10      = v10_p      ,       &
           ust     = ust_p       , xlat     = xlat_p     , xlong    = xlon_p     ,       &
           tskin   = tsk_p       , pblh     = hpbl_p     , t2m      = t2m_p      ,       &
           p8w     = pres2_hyd_p , dz8w     = dz_p       , z_at_w   = zgrid      ,       &
           p_phy   = pres_hyd_p  , t_phy    = t_p        , u_phy    = u_p        ,       &
           v_phy   = v_p         , qv       = qv_p       , vvel     = w_p        ,       &
           pi_phy  = pi_p        , rho_phy  = rho_p      , kpbl     = kpbl_p     ,       &
           nsoil   = num_soils   , smois    = smois_p    , tslb     = tslb_p     ,       &
           ivgtyp  = ivgtyp_p    , isltyp   = isltyp_p   , nlcat    = num_landcat,       &
           swdown  = swdown_p    , z0       = z0_p       , snowh    = snowh_p    ,       &
           julian  = curr_julday , rmol     = rmol_p     , raincv   = raincv_p   ,       &
           rainncv = rainncv_p   , dpt2m    = dpt2m_p    , znt      = znt_p      ,       &                              
           mavail  = mavail_p    , g        = gravity    ,                               &
           cp      = cp          , rd       = R_d        , gmt      = gmt        ,       &
           ids = ids , ide = ide , jds = jds , jde = jde , kds = kds , kde = kde ,       &
           ims = ims , ime = ime , jms = jms , jme = jme , kms = kds , kme = kme ,       &
           its = its , ite = ite , jts = jts , jte = jte , kts = kts , kte = kte         &
                            )
 call mpas_timer_stop('mpas_smoke')

!copy local arrays to MPAS grid:
 call smoke_to_MPAS(configs,diag_physics,tend_physics,its,ite)


 end subroutine driver_smoke

!=================================================================================================================
 end module mpas_atmphys_driver_smoke
!=================================================================================================================
