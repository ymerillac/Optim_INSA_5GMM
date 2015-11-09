! Copyright (C) 2010-2012  Daniel Filkovic

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! %                                                      %
! %            APAME - Aircraft Panel Method             %
! %______________________________________________________%
! %                                                      %
! %              3D potential flow solver                %
! %                                                      %
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! This file is part of APAME.

! APAME is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.

! APAME is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.

! You should have received a copy of the GNU General Public License
! along with APAME.  If not, see <http://www.gnu.org/licenses/>.

! file apame.f90

program apame

use module_input
use module_grid
use module_influence
use module_rhs
use module_velocity
use module_pressure

implicit none

! DECLARATIONS =================================================================
integer                                   :: interactive                    ,&
                                           & info_solver                    ,&
                                           & slash_position                 ,&
                                           & alloc_stat
real                                      :: pi                             ,&
                                           & start_time
integer,      allocatable, dimension(:,:) :: ipiv
real,         allocatable, dimension(:)   :: speed_x,speed_y,speed_z
character(len=100)                        :: file_name, job_name
character(len=3)                          :: version, subversion
character(len=6)                          :: build_date
logical                                   :: file_exists

parameter (pi=4.*atan(1.))
parameter (version="3")         ! program version
parameter (subversion="1")      ! program subversion
parameter (build_date="140915") ! build date
! ==============================================================================

! READING APAME INPUT FILE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! getting command line argument
call get_command_argument(1,file_name)
! in case no input file is given as argument, APAME will ask for one and enter interactive mode
if (file_name .eq. "") then
    ! entering interactive mode
    interactive=1
    
    ! printing title
    call func_title( version                                                ,&
                   & subversion                                             ,&
                   & build_date                                             )
    
    ! requesting user for file name
    read *, file_name
else
    ! entering non-interactive mode CHANGE THIS VALUE TO "0" WHEN RELEASING
    interactive=0
endif

call func_tic(start_time)
call func_new_line(interactive)

! check if input file exists
inquire(file=trim(file_name)//".inp",exist=file_exists)

! create (replace) log file
open(unit=2,file=trim(file_name)//".log",status="replace")

! open input file for usage
if (file_exists .neqv. .true.) then
    call func_message( interactive                                          ,&
                     & "    ERROR: Input file "//trim(file_name)//".inp not found")
    goto 999
endif

! print date and time
call subr_write_time(interactive)
call func_new_line(interactive)

! reading input file
call input(interactive, file_name, version, subversion)
if (input_err .ne. 0) then
    goto 999
endif

call func_message( interactive                                              ,&
                 & "Done!")
call func_new_line(interactive)
call func_new_line(interactive)

! calculating airspeed components
allocate ( speed_x(case_num)                                                ,&
         & speed_y(case_num)                                                ,&
         & speed_z(case_num)                                                )
call subr_speeds( case_num                                                  ,&
                & speed                                                     ,&
                & alfa,beta                                                 ,&
                & speed_x,speed_y,speed_z                                   )

! print jobname
job_name=trim(file_name)
! find slash or backslash position if file was given with path
slash_position=scan(job_name, '\/', back=.true.)
if (slash_position .eq. 0) then
    call func_message( interactive                                          ,&
                     & "    Job name: "//trim(job_name))
else
    call func_message( interactive                                          ,&
                     & "    Job name: "//job_name(slash_position+1:len_trim(job_name)))
endif

call func_new_line(interactive)
call func_new_line(interactive)

! print number of panels
call subr_write_ptot(interactive, panel_num, panel_num_no_wake)

! print keywords
call subr_write_key( interactive                                            ,&
                   & speed                                                  ,&
                   & ro                                                     ,&
                   & p_ref                                                  ,&
                   & mach                                                   ,&
                   & wing_span                                              ,&
                   & mac                                                    ,&
                   & wing_surf                                              ,&
                   & origin                                                 ,&
                   & method                                                 ,&
                   & error                                                  ,&
                   & coplanarity_angle                                      ,&
                   & farfield                                               ,&
                   & collcalc                                               ,&
                   & velorder                                               )

! print cases (angles of attack and sideslip angles)
call subr_write_case(interactive, case_num, alfa, beta)
call func_new_line(interactive)

! CALCULATING GRID INFORMATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

call func_message( interactive                                              ,&
                 & "Calculating grid information...")

call grid( interactive                                                      ,&
         & collcalc                                                         ,&
         & panel_num                                                        ,&
         & node_num                                                         ,&
         & panel_type                                                       ,&
         & node1,node2,node3,node4                                          ,&
         & farfield                                                         ,&
         & error                                                            ,&
         & colldist                                                         ,&
         & x_node,y_node,z_node                                             )
if (grid_err .ne. 0) then
    goto 999
endif

call func_message( interactive                                              ,&
                 & "Done!")
call func_new_line(interactive)
call func_new_line(interactive)

call func_toc( interactive, start_time                                      ,&
             & "    Time for preprocessing.............................:")
call func_new_line(interactive)


! CALCULATING INFLUENCE COEFFICIENTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

call func_tic(start_time)
call func_message( interactive                                              ,&
                 & "Calculating influence coefficients...")

call influence( method                                                      ,&
              & panel_num                                                   ,&
              & panel_num_no_wake                                           ,&
              & interactive                                                 ,&
              & panel_type                                                  ,&
              & elem1,elem2                                                 ,&
              & FF                                                          ,&
              & S                                                           ,&
              & error                                                       ,&
              & colx,coly,colz                                              ,&
              & cx,cy,cz                                                    ,&
              & l1,l2,l3                                                    ,&
              & p1,p2,p3                                                    ,&
              & n1,n2,n3                                                    ,&
              & x1,x2,x3,x4                                                 ,&
              & y1,y2,y3,y4                                                 ,&
              & d1,d2,d3,d4                                                 )
if (influence_err .ne. 0) then
    goto 999
endif

deallocate( d1,d2,d3,d4                                                     ,&
          & x1,x2,x3,x4                                                     ,&
          & y1,y2,y3,y4)

call func_message( interactive                                              ,&
                 & "Done!")
call func_new_line(interactive)
call func_new_line(interactive)

call func_toc( interactive, start_time                                      ,&
             & "    Time for calculating influence coefficients........:")
call func_new_line(interactive)

! CALCULATING RIGHT-HAND SIDE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

call func_tic(start_time)
call func_message( interactive                                              ,&
                 & "Calculating right-hand side...")

! allocate b to minimum size for constant doublet method
if (method .eq. 1) allocate (b(1,1))

call rhs_calc( b                                                            ,&
             & method                                                       ,&
             & panel_num                                                    ,&
             & panel_num_no_wake                                            ,&
             & interactive                                                  ,&
             & case_num                                                     ,&
             & speed_x,speed_y,speed_z                                      ,&
             & n1,n2,n3                                                     ,&
             & cx,cy,cz                                                     )
if (rhs_err .ne. 0) then
    goto 999
endif

if (method .eq. 0) then
    deallocate(b)
endif

call func_message( interactive                                              ,&
                 & "Done!")
call func_new_line(interactive)
call func_new_line(interactive)

call func_toc( interactive ,start_time                                      ,&
             & "    Time for calculating right-hand side...............:")
call func_new_line(interactive)

! SOLVING SYSTEM OF EQUATIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

call func_tic(start_time)
call func_message( interactive                                              ,&
                 & "Solving system of equations...")

allocate ( ipiv(panel_num_no_wake, panel_num_no_wake)                       ,&
         & stat=alloc_stat                                                  )
if (alloc_stat .ne. 0) then
    call func_message( interactive                                          ,&
                     & "    ERROR: Not enough memory for allocating permutation matrix &
                     & (main system of equations)")
    goto 999
endif

! solving system of equations
call sgesv( panel_num_no_wake                                               ,&
          & case_num                                                        ,&
          & a                                                               ,&
          & panel_num_no_wake                                               ,&
          & ipiv                                                            ,&
          & rhs                                                             ,&
          & panel_num_no_wake                                               ,&
          & info_solver                                                     )

deallocate(a, ipiv)

call func_message( interactive                                              ,&
                 & "Done!")
call func_new_line(interactive)

if (info_solver .lt. 0) then
    call func_message( interactive                                          ,&
                     & "    ERROR: Illegal value detected in system solver")
    goto 999
elseif (info_solver .gt. 0) then
    call func_message( interactive                                          ,&
                     & "    ERROR: Factorization complete with singular solution")
    goto 999
endif

call func_new_line(interactive)
call func_toc( interactive, start_time                                      ,&
             & "    Time for solving system of equations...............:")
call func_new_line(interactive)

! POSTPROCESSING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

call func_tic(start_time)
call func_message( interactive                                              ,&
                 & "Postprocessing...")

! calculating induced velocities
call velocities( rhs                                                        ,&
               & node_num                                                   ,&
               & panel_num                                                  ,&
               & panel_num_no_wake                                          ,&
               & case_num                                                   ,&
               & velorder                                                   ,&
               & method                                                     ,&
               & coplanarity_angle                                          ,&
               & panel_type                                                 ,&
               & interactive                                                ,&
               & cx,cy,cz                                                   ,&
               & x_node,y_node,z_node                                       ,&
               & (/l1,l2,l3/)                                               ,&
               & (/p1,p2,p3/)                                               ,&
               & (/n1,n2,n3/)                                               ,&
               & speed                                                      ,&
               & speed_x,speed_y,speed_z                                    ,&
               & node1,node2,node3,node4                                    ,&
               & elem1,elem2,elem3,elem4                                    )
if (velo_err .ne. 0) then
    goto 999
endif

! calculating pressure coefficients
call pressure( interactive                                                  ,&
             & panel_num                                                    ,&
             & panel_num_no_wake                                            ,&
             & panel_type                                                   ,&
             & case_num                                                     ,&
             & alfa,beta                                                    ,&
             & mac                                                          ,&
             & wing_span                                                    ,&
             & wing_surf                                                    ,&
             & origin                                                       ,&
             & speed                                                        ,&
             & mach                                                         ,&
             & ro                                                           ,&
             & p_ref                                                        ,&
             & S                                                            ,&
             & v                                                            ,&
             & n1,n2,n3                                                     ,&
             & cx,cy,cz                                                     )
if (press_err .ne. 0) then
    goto 999
endif

call func_message( interactive                                              ,&
                 & "Done!")
call func_new_line(interactive)
call func_new_line(interactive)
call func_toc( interactive, start_time                                      ,&
             & "    Time for postprocessing............................:")
call func_new_line(interactive)

! printing coefficients
call subr_write_coef( interactive                                           ,&
                    & case_num                                              ,&
                    & coef_x,coef_y,coef_z                                  ,&
                    & coef_l,coef_m,coef_n                                  ,&
                    & coef_drag,coef_side,coef_lift                         )

call func_new_line(interactive)
call func_new_line(interactive)

! WRITING APAME RESULTS FILE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (res_req .eq. 1) then

    ! allocate sigma to minimum size for constant doublet method
    if (method .eq. 1) allocate (sigma(1,1))

    call func_tic(start_time)
    call func_message( interactive                                          ,&
                     & "Writing results file...")
    
    call subr_results( interactive                                          ,&
                     & version                                              ,&
                     & subversion                                           ,&
                     & build_date                                           ,&
                     & file_name                                            ,&
                     & speed                                                ,&
                     & ro                                                   ,&
                     & p_ref                                                ,&
                     & mach                                                 ,&
                     & alfa                                                 ,&
                     & beta                                                 ,&
                     & wing_span                                            ,&
                     & mac                                                  ,&
                     & wing_surf                                            ,&
                     & origin                                               ,&
                     & method                                               ,&
                     & error                                                ,&
                     & coplanarity_angle                                    ,&
                     & farfield                                             ,&
                     & collcalc                                             ,&
                     & velorder                                             ,&
                     & case_num                                             ,&
                     & node_num                                             ,&
                     & panel_num                                            ,&
                     & panel_num_no_wake                                    ,&
                     & panel_type                                           ,&
                     & x_node,y_node,z_node                                 ,&
                     & node1,node2,node3,node4                              ,&
                     & elem1,elem2,elem3,elem4                              ,&
                     & cx,cy,cz                                             ,&
                     & cp                                                   ,&
                     & v                                                    ,&
                     & vx,vy,vz                                             ,&
                     & p_dyna                                               ,&
                     & p_mano                                               ,&
                     & p_stat                                               ,&
                     & rhs                                                  ,&
                     & sigma                                                ,&
                     & S                                                    ,&
                     & FF                                                   ,&
                     & n1,n2,n3                                             ,&
                     & l1,l2,l3                                             ,&
                     & p1,p2,p3                                             ,&
                     & coef_x,coef_y,coef_z                                 ,&
                     & coef_l,coef_m,coef_n                                 ,&
                     & coef_drag,coef_side,coef_lift                        ,&
                     & Fx,Fy,Fz                                             ,&
                     & Fl,Fm,Fn                                             ,&
                     & Fdrag,Fside,Flift                                    ,&
                     & results                                              )
    
    call func_message( interactive                                          ,&
                     & "Done!")
    call func_new_line(interactive)
    call func_new_line(interactive)
    call func_toc( interactive, start_time                                  ,&
                 & "    Time for writing results file......................:")
    call func_new_line(interactive)
endif

call func_message( interactive                                              ,&
                 & "Job complete!")
call func_new_line(interactive)

999 continue

call func_new_line(interactive)

close (2)

end program apame
