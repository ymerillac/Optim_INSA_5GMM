! Copyright (C) 2011  Daniel Filkovic

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

! file func_vect_inprod.f90

! This function calculates scalar product C of two vectors A and B.

subroutine func_vect_inprod(C, A, B)

implicit none

! INTENTS IN ===================================================================
real, dimension(3), intent(in)  :: A,B
! INTENTS OUT ==================================================================
real,               intent(out) :: C
! ==============================================================================

C = A(1)*B(1) + A(2)*B(2) + A(3)*B(3)

return
end subroutine func_vect_inprod
