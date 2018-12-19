
 BEGIN_PROVIDER [double precision,Energy_c_md_LDA_mu_of_r_grille_sphe_of_r,(N_states,n_points_radial_grid_spherical)]
&BEGIN_PROVIDER [double precision,Energy_c_md_LDA_mu_of_r_grille_sphe_of_r_integrated,(N_states)]
&BEGIN_PROVIDER [double precision, rho_r,(N_states,n_points_radial_grid_spherical)]
&BEGIN_PROVIDER [double precision,mu_r,(N_states,n_points_radial_grid_spherical)]
&BEGIN_PROVIDER [double precision,rhomu_r,(N_states,n_points_radial_grid_spherical)]
&BEGIN_PROVIDER [double precision,rhomu_r_integrated, (N_states)]
&BEGIN_PROVIDER [double precision,Energy_c_md_LDA_mu_of_r_integrated_core,(N_states)]
&BEGIN_PROVIDER [double precision,Energy_c_md_LDA_mu_of_r_integrated_valence,(N_states)]
 implicit none
 integer :: i,istate
 logical :: dospin
 double precision :: cpu0,cpu1
 call cpu_time(cpu0)
!!!!!!spherical integrale!!!!!
 rhomu_r_integrated = 0.d0
 Energy_c_md_LDA_mu_of_r_integrated_core(:) = 0.d0
 Energy_c_md_LDA_mu_of_r_integrated_valence(:) = 0.d0
 Energy_c_md_LDA_mu_of_r_grille_sphe_of_r_integrated = 0.d0 
 istate = 1
 Energy_c_md_LDA_mu_of_r_grille_sphe_of_r = 0.d0 
 rho_r = 0.d0
 mu_r = 0.d0
 rhomu_r = 0.d0
 ! initialization for OPENMP
 double precision :: ec_int,rho_r_int,mu_r_int,rhomu_r_int
 i = 1
 call give_all_spherical_averages(i,istate,ec_int,rho_r_int,mu_r_int,rhomu_r_int)

 !$OMP PARALLEL DO &
 !$OMP DEFAULT (NONE)  &
 !$OMP PRIVATE (i,ec_int,rho_r_int,rhomu_r_int,mu_r_int) & 
 !$OMP SHARED( n_points_radial_grid_spherical,istate,Energy_c_md_LDA_mu_of_r_grille_sphe_of_r,rho_r,rhomu_r,mu_r)
 do i = 1, n_points_radial_grid_spherical
  call give_all_spherical_averages(i,istate,ec_int,rho_r_int,mu_r_int,rhomu_r_int)
  Energy_c_md_LDA_mu_of_r_grille_sphe_of_r(istate,i) = ec_int
  rho_r(istate,i) = rho_r_int
  mu_r(istate,i) =  mu_r_int
  rhomu_r(istate,i) =  rhomu_r_int
  mu_r(istate,i) = mu_r_int
 enddo
 !$OMP END PARALLEL DO

 do i = 1, n_points_radial_grid_spherical
   Energy_c_md_LDA_mu_of_r_grille_sphe_of_r_integrated(istate) += Energy_c_md_LDA_mu_of_r_grille_sphe_of_r(istate,i) * pas_grid_sphe
   rhomu_r_integrated(istate) += rhomu_r(istate,i) * pas_grid_sphe
   if(list_r(i) < r_core_sphe) then
     Energy_c_md_LDA_mu_of_r_integrated_core(istate) += Energy_c_md_LDA_mu_of_r_grille_sphe_of_r(istate,i) * pas_grid_sphe   
   Else
     Energy_c_md_LDA_mu_of_r_integrated_valence(istate) += Energy_c_md_LDA_mu_of_r_grille_sphe_of_r(istate,i) * pas_grid_sphe 
   endif
 enddo
 call cpu_time(cpu1)
 print*,'Time for Energy_c_md_LDA_mu_of_r_integrated_core :',cpu1-cpu0


END_PROVIDER


subroutine give_all_spherical_averages(i,istate,ec_int,rho_r_int,mu_r_int,rhomu_r_int)
 implicit none
 integer, intent(in) :: i,istate
 double precision, intent(out) :: ec_int,rho_r_int,mu_r_int,rhomu_r_int
 integer :: j
 double precision :: r(3),aos_array(ao_num),two_dm,two_dm_laplacian,total_dm,two_dm_HF,two_dm_laplacian_HF,total_dm_HF
 double precision :: local_potential,two_body_dm,integral_f,mu,mu0,alpha,mu_integral
 double precision :: pi,dpi,rho_a(istate),rho_b(istate),ec(istate)
 logical :: dospin
 pi = dacos(-1.d0)
 dpi = 1.5d0 * dsqrt(pi)
 dospin = .True. ! JT dospin have to be set to true for open shell
 ec_int = 0.d0
 rho_r_int = 0.d0
 mu_r_int = 0.d0
 rhomu_r_int = 0.d0
   do j = 1, n_points_integration_angular
    r(1) = spherical_grid_of_r(1,j,i)
    r(2) = spherical_grid_of_r(2,j,i)
    r(3) = spherical_grid_of_r(3,j,i)

    call dm_dft_alpha_beta_and_all_aos_at_r(r,rho_a,rho_b,aos_array)
    if(mu_of_r_potential.EQ."cusp_condition")then
     call spherical_averaged_two_dm_at_second_order(r,0.d0,istate,two_dm,two_dm_laplacian,total_dm)
     call spherical_averaged_two_dm_HF_at_second_order(r,0.d0,istate,two_dm_HF,two_dm_laplacian_HF,total_dm_HF)
     two_dm = max(two_dm,1.d-15)
     two_dm_HF = max(two_dm_HF,1.d-15)
     mu0 =  dpi * (two_dm_laplacian / two_dm - two_dm_laplacian_HF / two_dm_HF)
     mu0 = max(mu0,1.d-15)
!    mu = mu0
     alpha = 1.d0 + 2.d0/(dsqrt(dacos(-1.d0)) * mu0)
     double precision :: alpha_bis, beta, delta
     alpha_bis = (two_dm_laplacian / two_dm - two_dm_laplacian_HF / two_dm_HF)
     alpha_bis = max(alpha_bis,1.d-15)
     beta = 2.d0/(3.d0*dacos(-1.d0))
     delta = 2.d0/dsqrt(dacos(-1.d0))
     mu = 1.d0/(2.d0*beta)*(alpha_bis + dsqrt(alpha_bis*alpha_bis + 4.d0 * alpha_bis * beta * delta))
    else if(mu_of_r_potential.EQ."integral_hf")then
 !   call integral_of_f_12_on_hf(r,integral_f)
 !   mu = mu_integral(integral_f,r)
    else if(mu_of_r_potential.EQ."hf_coallescence")then
     call f_HF_aa(r,r,local_potential)
     mu =  local_potential * dsqrt(dacos(-1.d0)) * 0.5d0
    else if(mu_of_r_potential.EQ."psi_coallescence")then
     call expectation_value_in_real_space(r,r,local_potential,two_body_dm)
     mu =  local_potential * dsqrt(dacos(-1.d0)) * 0.5d0
    else 
      print*,'you requested the following mu_of_r_potential'
      print*,mu_of_r_potential
      print*,'which does not correspond to any of the options for such keyword'
      stop
    endif
    if(mu.lt.0.d0)then
     print*,r
     print*,mu
    endif
!!!!!!!!!! CORRELATION PART
    call ESRC_MD_LDAERF (mu,rho_a(istate),rho_b(istate),dospin,ec(istate))
    ec_int += weights_rad_of_r(j,i) * ec(istate) 
    rho_r_int += (rho_a(istate)+rho_b(istate)) * weights_rad_of_r(j,i)
    mu_r_int +=  mu * weights_rad_of_r(j,i) 
    rhomu_r_int +=  (rho_a(istate)+rho_b(istate)) * mu * weights_rad_of_r(j,i)       
    rhomu_r_int = rhomu_r_int / dble(elec_num)
  enddo
    if(list_r(i)*list_r(i) * 4.d0 * pi.gt.1.d-10)then
     mu_r_int = mu_r_int  / (list_r(i)*list_r(i) * 4.d0 * pi)
    else
     mu_r_int = 0.d0
    endif

end

!   Energy_c_md_LDA_mu_of_r_grille_sphe_of_r(istate,i) += weights_rad_of_r(j,i) * ec(istate) 
!   rho_r(istate,i) += (rho_a(istate)+rho_b(istate)) * weights_rad_of_r(j,i)
!   mu_r(istate,i) +=  mu * weights_rad_of_r(j,i) 
!   rhomu_r(istate,i) +=  (rho_a(istate)+rho_b(istate)) * mu * weights_rad_of_r(j,i)       
!   rhomu_r(istate,i) = rhomu_r(istate,i) / dble(elec_num)
!   mu_r(istate,i) = mu_r(istate,i) /  max(list_r(i)*list_r(i) * 4.d0 * pi,1.d-10)
! enddo
