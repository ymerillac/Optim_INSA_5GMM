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

! file func_message.f90

! This function prints text string on screen and/or log file depending on
! integer interactive.

subroutine func_message(interactive, string)

implicit none

! INTENTS IN ===================================================================
integer,          intent(in) :: interactive
character(len=*), intent(in) :: string
! ==============================================================================

if (interactive .eq. 1) then
    write (6,100,advance="no") string
    flush (6)
    write (2,100,advance="no") string
    flush (2)
else
    write (2,100,advance="no") string
    flush (2)
endif

! FORMATS ======================================================================
100 format(1x,a)
! ==============================================================================

return
end subroutine func_message
