!> Contains routines necessary to initialize the Database client
module MOM_database_comms

use MOM_cpu_clock,        only : cpu_clock_id, cpu_clock_begin, cpu_clock_end, CLOCK_ROUTINE
use MOM_error_handler,    only : MOM_error, FATAL, WARNING, NOTE, MOM_mesg, is_root_pe
use MOM_file_parser,      only : read_param, get_param, log_version, param_file_type
use smartredis_client,    only : dbclient_type => client_type

implicit none; private

!> Control structure to store Database client related parameters and objects
type, public :: dbcomms_CS_type
  type(dbclient_type) :: client !< The Database client itself
  logical           :: use_dbclient !< If True, use Database within MOM6
  logical           :: colocated !< If True, the orchestrator was setup in 'co-located' mode
  logical           :: cluster   !< If True, the orchestrator has three shards or more
  integer           :: colocated_stride !< Sets which ranks will load the model from the file
                                        !! e.g. mod(rank,colocated_stride) == 0
end type dbcomms_CS_type

public :: dbclient_type
public :: database_comms_init

contains

subroutine database_comms_init(param_file, CS, client_in)
  type(param_file_type),       intent(in   ) :: param_file !< Parameter file structure
  type(dbcomms_CS_type),    intent(inout) :: CS         !< Control structure for Database
  type(dbclient_type), optional, intent(in   ) :: client_in !< If present, use a previously initialized
                                                          !! Database client

  character(len=40) :: mdl = "MOM_DBCLIENT"
  integer :: id_client_init
  integer :: return_code
  call get_param(param_file, mdl, "USE_DBCLIENT",  CS%use_dbclient, &
                 "If true, use the data client to connect"//&
                 "with the Database database", default=.false.)

  if (present(client_in)) then ! The driver (e.g. the NUOPC cap) has already initialized the client

    CS%client = client_in

    if (.not. CS%client%isinitialized() .and. CS%use_dbclient) then
      call MOM_error(FATAL, &
      "If using a Database client not initialized within MOM, client%initialize must have already been invoked."//&
      " Check that the client has been initialized in the driver before the call to initialize_MOM")
    endif

  elseif (CS%use_dbclient) then ! The client will be initialized within MOM

    call get_param(param_file, mdl, "DBCLIENT_COLOCATED",  CS%colocated, &
                   "If true, the Database database is colocated on the simulation nodes.",&
                   default=.false.)
    if (CS%colocated) then
      CS%cluster = .false.
      call get_param(param_file, mdl, "DBCLIENT_COLOCATED_STRIDE",  CS%colocated_stride, &
                     "If true, the Database database is colocated on the simulation nodes.",&
                     default=0)
    else
      call get_param(param_file, mdl, "DBCLIENT_CLUSTER",  CS%cluster, &
                     "If true, the Database database is distributed over multiple nodes.",&
                     default=.true.)
    endif
    id_client_init = cpu_clock_id('(DBCLIENT client init)', grain=CLOCK_ROUTINE)
    call MOM_error(NOTE,"Database Client Initializing")
    call cpu_clock_begin(id_client_init)
    return_code = CS%client%initialize(CS%cluster)
    if (CS%client%SR_error_parser(return_code)) then
      call MOM_error(FATAL, "Database client failed to initialize")
    endif
    call MOM_error(NOTE,"Database Client Initialized")
    call cpu_clock_end(id_client_init)

  endif
end subroutine database_comms_init

end module MOM_database_comms

