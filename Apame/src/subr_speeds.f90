! Copyright (C) 2010 Daniel Filkovic

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

! file subr_speeds.f90

! This function calculates free stream velocity components (speed_x,speed_y,speed_z)
! with given angles of attack (alfa), sideslip angles (beta) and total velocity
! (speed) for each case.

subroutine subr_speeds( case_num                                            ,&
                      & speed                                               ,&
                      & alfa,beta                                           ,&
                      & speed_x,speed_y,speed_z                             )

implicit none

! INTENTS IN ===================================================================
integer,                      intent(in)  :: case_num
real,                         intent(in)  :: speed
real,    dimension(case_num), intent(in)  :: alfa
real,    dimension(case_num), intent(in)  :: beta
! INTENTS OUT ==================================================================
real,    dimension(case_num), intent(out) :: speed_x,speed_y,speed_z
! PRIVATE ======================================================================
integer                                   :: i
! ==============================================================================

do i=1,case_num
    speed_x(i) = speed*cos(alfa(i))*cos(beta(i))
    speed_y(i) = speed*sin(beta(i))
    speed_z(i) = speed*cos(beta(i))*sin(alfa(i))
enddo

return
end subroutine subr_speeds
