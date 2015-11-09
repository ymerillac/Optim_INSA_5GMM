! Copyright (C) 2010-2011  Daniel Filkovic

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

! file func_tic.f90

! This function returns current time (start_time) in seconds.

subroutine func_tic(start_time)

implicit none

! INTENTS OUT ==================================================================
real,          intent(out) :: start_time
! PRIVATE ======================================================================
integer                    :: hour,minute
real                       :: second 
character(10)              :: time
! ==============================================================================

! getting current time
call date_and_time(TIME=time)

! reading hour minute and second from "time"
read (time,100) hour
read (time,101) minute
read (time,102) second

! calculating time in seconds
start_time = hour*3600 + minute*60 + second

! FORMATS ======================================================================
100 format(i2)
101 format(2x,i2)
102 format(4x,F6.3)
! ==============================================================================

return
end subroutine func_tic
