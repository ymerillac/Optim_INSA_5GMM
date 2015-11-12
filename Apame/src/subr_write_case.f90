! Copyright (C) 2010  Daniel Filkovic

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

! file subr_write_case.f90

! This subroutine prints on screen and/or log file requested angles of attack
! (alfa) and sideslip angles (beta) in degrees depending on integer interactive.

subroutine subr_write_case(interactive, case_num, alfa, beta)

implicit none

! INTENTS IN ===================================================================
integer,                      intent(in) :: interactive, case_num
real,    dimension(case_num), intent(in) :: alfa
real,    dimension(case_num), intent(in) :: beta
! PRIVATE ======================================================================
integer                                  :: i
real                                     :: pi
real,    dimension(case_num)             :: alfa_deg
real,    dimension(case_num)             :: beta_deg
parameter (pi=4.*atan(1.))
! ==============================================================================

! converting angles to degrees
alfa_deg = alfa/pi*180
beta_deg = beta/pi*180

! printing angles of attack
call func_message( interactive                                              ,&
                 & "    Angles of interest [deg]:")
call func_new_line(interactive)
call func_message( interactive                                              ,&
                 & "        Angles of attack:")
if (interactive .eq. 1) then
    do i=1,case_num
        write (6,100,advance="no") alfa_deg(i)
        write (2,100,advance="no") alfa_deg(i)
    enddo
else
    do i=1,case_num
        write (2,100,advance="no") alfa_deg(i)
    enddo
endif

! printing sideslip angles
call func_new_line(interactive)
call func_message( interactive                                              ,&
                 & "        Sideslip angles: ")
if (interactive .eq. 1) then
    do i=1,case_num
        write (6,100,advance="no") beta_deg(i)
        write (2,100,advance="no") beta_deg(i)
    enddo
else
    do i=1,case_num
        write (2,100,advance="no") beta_deg(i)
    enddo
endif
call func_new_line(interactive)

! FORMATS ======================================================================
100 format(1x,f5.1)
! ==============================================================================

return
end subroutine subr_write_case
