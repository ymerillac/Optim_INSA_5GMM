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

! file subr_write_key.f90

! This subroutine prints on screen and/or log file keywords and their values
! from input file depending on integer interactive.

subroutine subr_write_key( interactive                                      ,&
                         & speed                                            ,&
                         & ro                                               ,&
                         & p_ref                                            ,&
                         & mach                                             ,&
                         & wing_span                                        ,&
                         & mac                                              ,&
                         & wing_surf                                        ,&
                         & origin                                           ,&
                         & method                                           ,&
                         & error                                            ,&
                         & coplanarity_angle                                ,&
                         & farfield                                         ,&
                         & collcalc                                         ,&
                         & velorder                                         )

implicit none

! INTENTS IN ===================================================================
integer,               intent(in) :: interactive                            ,&
                                   & collcalc                               ,&
                                   & velorder                               ,&
                                   & method
real,                  intent(in) :: speed                                  ,&
                                   & ro                                     ,&
                                   & p_ref                                  ,&
                                   & mach                                   ,&
                                   & wing_span                              ,&
                                   & mac                                    ,&
                                   & wing_surf                              ,&
                                   & error                                  ,&
                                   & coplanarity_angle                      ,&
                                   & farfield
real,    dimension(3), intent(in) :: origin
! PRIVATE ======================================================================
real                              :: pi

parameter (pi=4.*atan(1.))
! ==============================================================================

call func_message( interactive                                              ,&
                 & "    KEYWORD values:")
call func_new_line(interactive)

if (interactive .eq. 1) then
    write (6,100) "Airspeed:                 ",speed
    write (6,100) "Density:                  ",ro
    write (6,100) "Atmospheric pressure:     ",p_ref
    write (6,100) "Mach number:              ",mach
    write (6,100) "Wing span:                ",wing_span
    write (6,100) "Mean aerodynamic chord:   ",mac
    write (6,100) "Wing surface:             ",wing_surf
    write (6,101) "Origin:                   ",origin
    if (method .eq. 0) then
        write (6,102) "Singularity method:       constant source/doublet"
    elseif (method .eq. 1) then
        write (6,102) "Singularity method:       constant doublet"
    endif
    write (6,100) "Error parameter:          ",error
    write (6,100) "Farfield parameter:       ",farfield
    if (collcalc .eq. 0) then
        write (6,102) "Collocation calculation:  approximative"
    elseif (collcalc .eq. 1) then
        write (6,102) "Collocation calculation:  accurate"
    endif
    if (velorder .eq. 0) then
        write (6,102) "Velocity interp. method:  nodal"
    elseif (velorder .eq. 1) then
        write (6,102) "Velocity interp. order:   linear"
        write (6,100) "Curvature corr. angle:    ",coplanarity_angle*180/pi
    elseif (velorder .eq. 2) then
        write (6,102) "Velocity interp. order:   quadratic"
        write (6,100) "Curvature corr. angle:    ",coplanarity_angle*180/pi
    endif
    write (2,100) "Airspeed:                 ",speed
    write (2,100) "Density:                  ",ro
    write (2,100) "Atmospheric pressure:     ",p_ref
    write (2,100) "Mach number:              ",mach
    write (2,100) "Wing span:                ",wing_span
    write (2,100) "Mean aerodynamic chord:   ",mac
    write (2,100) "Wing surface:             ",wing_surf
    write (2,101) "Origin:                   ",origin
    if (method .eq. 0) then
        write (2,102) "Singularity method:       constant source/doublet"
    elseif (method .eq. 1) then
        write (2,102) "Singularity method:       constant doublet"
    endif
    write (2,100) "Error parameter:          ",error
    write (2,100) "Farfield parameter:       ",farfield
    if (collcalc .eq. 0) then
        write (2,102) "Collocation calculation:  approximative"
    elseif (collcalc .eq. 1) then
        write (2,102) "Collocation calculation:  accurate"
    endif
    if (velorder .eq. 0) then
        write (2,102) "Velocity interp. method:  nodal"
    elseif (velorder .eq. 1) then
        write (2,102) "Velocity interp. order:   linear"
        write (2,100) "Curvature corr. angle:    ",coplanarity_angle*180/pi
    elseif (velorder .eq. 2) then
        write (2,102) "Velocity interp. order:   quadratic"
        write (2,100) "Curvature corr. angle:    ",coplanarity_angle*180/pi
    endif
else
    write (2,100) "Airspeed:                 ",speed
    write (2,100) "Density:                  ",ro
    write (2,100) "Atmospheric pressure:     ",p_ref
    write (2,100) "Mach number:              ",mach
    write (2,100) "Wing span:                ",wing_span
    write (2,100) "Mean aerodynamic chord:   ",mac
    write (2,100) "Wing surface:             ",wing_surf
    write (2,101) "Origin:                   ",origin
    if (method .eq. 0) then
        write (2,102) "Singularity method:       constant source/doublet"
    elseif (method .eq. 1) then
        write (2,102) "Singularity method:       constant doublet"
    endif
    write (2,100) "Error parameter:          ",error
    write (2,100) "Farfield parameter:       ",farfield
    if (collcalc .eq. 0) then
        write (2,102) "Collocation calculation:  approximative"
    elseif (collcalc .eq. 1) then
        write (2,102) "Collocation calculation:  accurate"
    endif
    if (velorder .eq. 0) then
        write (2,102) "Velocity interp. method:  nodal"
    elseif (velorder .eq. 1) then
        write (2,102) "Velocity interp. order:   linear"
        write (2,100) "Curvature corr. angle:    ",coplanarity_angle*180/pi
    elseif (velorder .eq. 2) then
        write (2,102) "Velocity interp. order:   quadratic"
        write (2,100) "Curvature corr. angle:    ",coplanarity_angle*180/pi
    endif
endif
call func_new_line(interactive)

! FORMATS ======================================================================
100 format(9x,a,e10.5)
101 format(9x,a,e10.5,1x,e10.5,1x,e10.5)
102 format(9x,a)
! ==============================================================================

return
end subroutine subr_write_key
