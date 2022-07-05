!> Contains routines that contain dummy routines for the
module MOM_MEKE_smartredis

use iso_c_binding,         only : c_float

use MOM_coms,              only : PE_here
use MOM_cpu_clock,         only : cpu_clock_id, cpu_clock_begin, cpu_clock_end, CLOCK_ROUTINE
use MOM_diag_mediator,     only : post_data, register_diag_field
use MOM_diag_mediator,     only : diag_ctrl, time_type
use MOM_domains,           only : pass_var, pass_vector
use MOM_error_handler,     only : MOM_error, FATAL, WARNING, is_root_pe
use MOM_file_parser,       only : read_param, get_param, log_version, param_file_type
use MOM_grid,              only : ocean_grid_type
use MOM_interface_heights, only : find_eta
use MOM_isopycnal_slopes,  only : calc_isoneutral_slopes
use MOM_smartredis,        only : client_type, smartredis_CS_type
use MOM_time_manager,      only : time_type_to_real
use MOM_variables,         only : thermo_var_ptrs
use MOM_unit_scaling,      only : unit_scale_type
use MOM_verticalGrid,      only : verticalGrid_type

implicit none; private

#include <MOM_memory.h>
#include <enum_fortran.inc>

public meke_smartredis_init, infer_meke

type, public :: meke_smartredis_CS_type; private

  ! Pointers to derived types from other modules
  type(client_type), pointer :: client => NULL() !< Pointer to the SmartRedis client
  type(diag_ctrl),   pointer :: diag => NULL() !< A type that regulates diagnostics output

  ! Constants for this module
  integer, parameter :: num_features = 4 !< How many features used to predict EKE
  integer, parameter :: mke_idx = 1   !< Index of mean kinetic energy in the feature array
  integer, parameter :: slope_z_idx = !< Index of vertically averaged isopycnal slope in the feature array
  integer, parameter :: rv_idx =      !< Index of surface relative vorticity in the feature array
  integer, parameter :: rd_dx_z_idx = !< Index of the radius of deformation over the grid size in the feature array

  ! Inferring EKE from ML
  logical :: online_analysis !< If true, post the EKE used in MOM6 at every timestep
  character(len=5) :: model_key  = 'mleke'  !< Key where the ML-model is stored
  character(len=7) :: key_suffix !< Suffix appended to every key sent to Redis
  real :: log_eke_max ! The maximum value of log(eke) considered physically reasonable

  ! Clock ids
  integer :: id_client_init   !< Clock id to time initialization of the client
  integer :: id_put_tensor    !< Clock id to time put_tensor routine
  integer :: id_run_model     !< Clock id to time running of the ML model
  integer :: id_unpack_tensor !< Clock id to time retrieval of EKE prediction

  ! Diagnostic ids
  integer :: id_mke     = -1 !< Diagnostic id for surface mean kinetic energy
  integer :: id_slope_z = -1 !< Diagnostic id for vertically averaged horizontal slope magnitude
  integer :: id_slope_x = -1 !< Diagnostic id for isopycnal slope in the x-direction
  integer :: id_slope_y = -1 !< Diagnostic id for isopycnal slope in the y-direction
  integer :: id_rv      = -1 !< Diagnostic id for surface relative vorticity

end type meke_smartredis_CS_type

contains

!> Initializer for the SmartRedis MEKE module that uses ML to predict eddy kinetic energy
subroutine meke_smartredis_init(diag, G, US, param_file, smartredis_CS, CS)
  type(diag_ctrl), target, intent(inout) :: diag       !< Diagnostics structure.
  type(ocean_grid_type),         intent(inout) :: G          !< The ocean's grid structure.
  type(unit_scale_type),         intent(in)    :: US         !< A dimensional unit scaling type
  type(param_file_type),         intent(in)    :: param_file !< Parameter file parser structure.
  type(smartredis_CS_type),      intent(in)    :: smartredis_CS !< SmartRedis client
  type(meke_smartredis_CS_type), intent(inout) :: CS         !< Control structure for this module

  character(len=200)  :: inputdir
  integer :: sr_return_code

  ! Store pointers in control structure
  CS%diag => diag
  CS%client => smartedis_CS%client

  write(CS%key_suffix, '(A,I6.6)') '_', PE_here()
  ! Put some basic information into the database
  sr_return_code = 0
  sr_return_code = CS%client%put_tensor("meta"//CS%key_suffix, &
    REAL([G%isd_global, G%idg_offset, G%jsd_global, G%jdg_offset]),[4]) + sr_return_code
  sr_return_code = CS%client%put_tensor("geolat"//CS%key_suffix, G%geoLatT, shape(G%geoLatT)) + sr_return_code
  sr_return_code = CS%client%put_tensor("geolon"//CS%key_suffix, G%geoLonT, shape(G%geoLonT)) + sr_return_code
  sr_return_code = CS%client%put_tensor("EKE_shape"//CS%key_suffix, shape(MEKE%MEKE), [2]) + sr_return_code

  if (sr_return_code /= SRNoError) call MOM_error(FATAL, "Putting metadata into the database failed")

  call read_param(param_file, "INPUTDIR", inputdir)
  inputdir = slasher(inputdir)

  call get_param(param_file, mdl, "BATCH_SIZE", batch_size, "Batch size to use for inference", default=1)
  call get_param(param_file, mdl, "EKE_BACKEND", backend, &
                 "The computational backend to use for EKE inference (CPU or GPU)", default="GPU")
  call get_param(param_file, mdl, "EKE_MODEL", model_filename, &
                 "Filename of the a saved pyTorch model to use", fail_if_missing = .true.)
  call get_param(param_file, mdl, "EKE_MAX", CS%eke_max, &
                 "Maximum value of EKE allowed when inferring EKE", default=2., scale=US%L_T_to_m_s**2)

  ! Set the machine learning model
  if (smartredis_CS%colocated) then
    if (modulo(PE_here(),smartredis_CS%colocated_stride) == 0) then
      sr_return_code = CS%client%set_model_from_file(CS%model_key, trim(CS%inputdir)//trim(model_filename), &
                                                  "TORCH", backend, batch_size=batch_size)
    endif
  else
    if (is_root_pe()) then
      sr_return_code = CS%client%set_model_from_file(CS%model_key, trim(CS%inputdir)//trim(model_filename), &
                                                  "TORCH", backend, batch_size=batch_size)
      if (CS%client%SR_error_parser(sr_return_code)) then
        call MOM_error(FATAL, "MEKE: set_model failed")
      endif
    endif
  endif
  if (sr_return_code /= SRNoError) call MOM_error(FATAL,"Error loading machine learning model into the database")

  call get_param(param_file, mdl, "ONLINE_ANALYSIS", CS%online_analysis, &
               "If true, post EKE used in MOM6 to the database for analysis", default=.true.)

  ! Set various clock ids
  CS%id_client_init   = cpu_clock_id('(SMARTREDIS client init)', grain=CLOCK_ROUTINE)
  CS%id_put_tensor    = cpu_clock_id('(SMARTREDIS put tensor)', grain=CLOCK_ROUTINE)
  CS%id_run_model     = cpu_clock_id('(SMARTREDIS run model)', grain=CLOCK_ROUTINE)
  CS%id_run_script    = cpu_clock_id('(SMARTREDIS run script)', grain=CLOCK_ROUTINE)
  CS%id_unpack_tensor = cpu_clock_id('(SMARTREDIS unpack tensor )', grain=CLOCK_ROUTINE)

  ! Diagnostics for SMARTREDIS
  CS%id_mke = register_diag_field('ocean_model', 'MEKE_MKE', diag%axesT1, Time, &
     'Surface mean (resolved) kinetic energy used in MEKE', 'm2 s-2', conversion=US%L_T_to_m_s**2)
  CS%id_slope_z= register_diag_field('ocean_model', 'MEKE_slope_z', diag%axesT1, Time, &
     'Vertically averaged isopyncal slope magnitude used in MEKE', 'm2 s-2', conversion=US%L_T_to_m_s**2)
  CS%id_slope_x= register_diag_field('ocean_model', 'MEKE_slope_x', diag%axesCui, Time, &
     'Isopycnal slope in the x-direction used in MEKE', 'm2 s-2', conversion=US%L_T_to_m_s**2)
  CS%id_slope_y= register_diag_field('ocean_model', 'MEKE_slope_y', diag%axesCvi, Time, &
     'Isopycnal slope in the y-direction used in MEKE', 'm2 s-2', conversion=US%L_T_to_m_s**2)
  CS%id_rv= register_diag_field('ocean_model', 'MEKE_RV', diag%axesT1, Time, &
     'Surface relative vorticity used in MEKE', 'm2 s-2', conversion=US%L_T_to_m_s**2)

end subroutine meke_smartredis_init

!> Use the SmartRedis client to call a machine learning to predict eddy kinetic energy
subroutine infer_meke(MEKE, u, v, tv, h, dt, G, GV, CS)
  real, dimension(SZI_(G),SZJ_(G)), intent(  out) :: MEKE !< Vertically averaged eddy kinetic energy [L2 T-2 ~> m2 s-2]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)), intent(inout) :: u  !< Zonal velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)), intent(inout) :: v  !< Meridional velocity [L T-1 ~> m s-1]
  type(thermo_var_ptrs),                     intent(in)    :: tv !< Type containing thermodynamic variables
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in)    :: h  !< Layer thickness [H ~> m or kg m-2].
  real,                                      intent(in)    :: dt !< Model(baroclinic) time-step [T ~> s].
  type(ocean_grid_type),                     intent(inout) :: G  !< Ocean grid
  type(verticalGrid_type),                   intent(in)    :: GV !< Ocean vertical grid structure
  type(meke_smartredis_CS_type),             intent(in)    :: CS !< Control structure for inferring MEKE using SmartRedis

  real, dimension(SZI_(G),SZJ_(G)) :: mke
  real, dimension(SZI_(G),SZJ_(G)) :: slope_z
  real, dimension(SZI_(G),SZJ_(G)) :: rv_z
  real, dimension(SZI_(G),SZJ_(G)) :: rd_dx_z

  real, dimension(SZIB_(G),SZJ_(G), SZK_(G)) :: h_u ! Thickness at u point
  real, dimension(SZI_(G),SZJB_(G), SZK_(G)) :: h_v ! Thickness at v point
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)+1) :: slope_x ! Isoneutral slope at U point
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)+1) :: slope_y ! Isoneutral slope at V point
  real, dimension(SZIB_(G),SZJ_(G)) :: slope_x_vert_avg ! Isoneutral slope at U point
  real, dimension(SZI_(G),SZJB_(G)) :: slope_y_vert_avg ! Isoneutral slope at V point
  real, dimension(SZI_(G), SZJ_(G), SZK_(G)+1) ::  e ! The interface heights relative to mean sea level [Z ~> m].

  character(len=255), dimension(1) :: model_out, model_in
  real(kind=c_float), dimension(SZI_(G)*SZJ_(G),CS%num_features) :: features_array
  real(kind=c_float), dimension(SZI_(G)*SZJ_(G)) :: MEKE_vec

  real :: slope_t, u_t, v_t ! u and v interpolated to thickness point
  real :: dvdx, dudy

  integer :: i, j, k, is, ie, js, je, Isq, Ieq, Jsq, Jeq, nz
  integer sr_return_code

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB

  ! Calculate various features for used to infer eddy kinetic energy
  ! Linear interpolation to estimate thickness at a velocity points
  do k=1,nz; do j=js-1,je+1; do i=is-1,ie+1
    h_u(I,j,k) = 0.5*(h(i,j,k)*G%mask2dT(i,j) + h(i+1,j,k)*G%mask2dT(i+1,j)) + GV%Angstrom_H
    h_v(i,J,k) = 0.5*(h(i,j,k)*G%mask2dT(i,j) + h(i,j+1,k)*G%mask2dT(i,j+1)) + GV%Angstrom_H
  enddo; enddo; enddo;
  call find_eta(h, tv, G, GV, US, e, halo_size=2)
  call calc_isoneutral_slopes(G, GV, US, h, e, tv, dt*1.e-7, slope_x, slope_y)
  call pass_vector(slope_x, slope_y, G%Domain)
  do j=js-1,je+1; do i=is-1,ie+1
    slope_x_vert_avg(I,j) = vertical_average_interface(slope_x(i,j,:), h_u(i,j,:), GV%H_subroundoff)
    slope_y_vert_avg(i,J) = vertical_average_interface(slope_y(i,j,:), h_v(i,j,:), GV%H_subroundoff)
  enddo; enddo
  slope_z(:,:) = 0.

  call pass_vector(slope_x_vert_avg, slope_y_vert_avg, G%Domain)
  call pass_vector(u, v, G%Domain)
  do j=js,je; do i=is,ie
    ! Calculate weights for interpolation from velocity points to h points
    sum_area = G%areaCu(I-1,j) + G%areaCu(I,j)
    if (sum_area>0.0) then
      Idenom = sqrt(0.5*G%IareaT(i,j) / sum_area)
      a_w = G%areaCu(I-1,j) * Idenom
      a_e = G%areaCu(I,j) * Idenom
    else
      a_w = 0.0 ; a_e = 0.0
    endif

    sum_area = G%areaCv(i,J-1) + G%areaCv(i,J)
    if (sum_area>0.0) then
      Idenom = sqrt(0.5*G%IareaT(i,j) / sum_area)
      a_s = G%areaCv(i,J-1) * Idenom
      a_n = G%areaCv(i,J) * Idenom
    else
      a_s = 0.0 ; a_n = 0.0
    endif

    ! Calculate mean kinetic energy
    u_t = a_e*u(I,j,1)+a_w(I-1,j,1)
    v_t = a_n*v(i,J,1)+a_s(i,J-1,1))
    mke(i,j) = 0.5*( u_t*u_t + v_t*v_t )

    ! Calculate the magnitude of the slope
    slope_t = slope_x_vert_avg(I,j)*a_e+slope_x_vert_avg(I-1,j)*a_w
    slope_z(i,j) = sqrt(slope_t*slope_t)
    slope_t = slope_y_vert_avg(i,J)*a_n+slope_y_vert_avg(i,J-1)*a_s
    slope_z(i,j) = 0.5*(slope_z(i,j) + sqrt(slope_t*slope_t))*G%mask2dT(i,j)
  enddo; enddo
  call pass_var(slope_z, G%Domain)
  call pass_var(MEKE%Rd_dx_h, G%Domain)

  ! Calculate relative vorticity
  do J=Jsq-1,Jeq+1 ; do I=Isq-1,Ieq+1
    dvdx = (v(i+1,J,1)*G%dyCv(i+1,J) - v(i,J,1)*G%dyCv(i,J))
    dudy = (u(I,j+1,1)*G%dxCu(I,j+1) - u(I,j,1)*G%dxCu(I,j))
    ! Assumed no slip
    rv_z(I,J) = (2.0-G%mask2dBu(I,J)) * (dvdx - dudy) * G%IareaBu(I,J)
  enddo; enddo

  ! Construct the feature array
  features_array(:,CS%mke_idx) = pack(mke,.true.)
  features_array(:,CS%slope_z_idx) = pack(slope_z,.true.)
  features_array(:,CS%rd_dx_z_idx) = pack(MEKE%Rd_dx_h,.true.)
  features_array(:,CS%rv_idx) = pack(rv_z,.true.)
  call cpu_clock_begin(CS%id_put_tensor)
  sr_return_code = client%put_tensor("features"//CS%key_suffix, features_array, shape(features_array))
  call cpu_clock_end(CS%id_put_tensor)

  ! Run the ML model to predict EKE and return the result
  model_out(1) = "EKE"//CS%key_suffix
  model_in(1) = "features"//CS%key_suffix
  call cpu_clock_begin(CS%id_run_model)
  sr_return_code = CS%client%run_model(CS%model_key, model_in, model_out)
  call cpu_clock_end(CS%id_run_model)
  if (client%SR_error_parser(sr_return_code)) then
    call MOM_error(FATAL, "MEKE: run_model failed")
  endif
  call cpu_clock_begin(CS%id_unpack_tensor)
  sr_return_code = client%unpack_tensor( model_out(1), MEKE_vec, shape(MEKE_vec) )
  call cpu_clock_end(CS%id_unpack_tensor)

  MEKE%MEKE = reshape(CS%MEKE_vec, shape(MEKE%MEKE))
  do j=js,je; do i=is,ie
    MEKE%MEKE(i,j) = MIN(MAX(exp(MEKE%MEKE(i,j)),0.),CS%log_eke_max)
  enddo; enddo
  call pass_var(MEKE%MEKE,G%Domain)

  if (CS%online_analysis) then
    write(time_suffix,"(F16.0)") time_type_to_real(Time)
    sr_return_code = client%put_tensor(trim("EKE_")//trim(adjustl(time_suffix))//CS%key_suffix, MEKE%MEKE, shape(MEKE%MEKE))
  endif

  if (CS%id_rv>0) call post_data(CS%id_rv, rv_z, CS%diag)
  if (CS%id_mke>0) call post_data(CS%id_mke, mke, CS%diag)
  if (CS%id_slope_z>0) call post_data(CS%id_slope_z, slope_z, CS%diag)
  if (CS%id_slope_x>0) call post_data(CS%id_slope_x, slope_x, CS%diag)
  if (CS%id_slope_y>0) call post_data(CS%id_slope_y, slope_y, CS%diag)
end subroutine infer_meke

!> Compute average of interface quantities weighted by the thickness of the surrounding layers
real function vertical_average_interface(h, w, h_min)

  real, dimension(:), intent(in) :: h  !< Layer Thicknesses
  real, dimension(:), intent(in) :: w !< Quantity to average
  real, intent(in) :: h_min !< The vanishingly small layer thickness

  real :: htot, inv_htot
  integer :: k, nk

  nk = size(uh)
  htot = h_min
  do k=2,nk
    htot = htot + (h(k-1)+h(k))
  enddo
  inv_htot = 1./htot

  vertical_average = 0.
  do K=2,nk
    vertical_average = vertical_average + (w(k)*(h(k-1)+h(k)))*inv_htot
  enddo
end function vertical_average

end module smartredis_meke