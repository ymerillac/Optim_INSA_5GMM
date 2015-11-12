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

! file module_pressure.f90

! This module calculates pressure distributions and sums forces

module module_pressure

! DECLARATIONS =================================================================
integer                           :: press_err
real, allocatable, dimension(:,:) :: cp                                     ,&
                                   & p_dyna                                 ,&
                                   & p_mano                                 ,&
                                   & p_stat
real, allocatable, dimension(:)   :: Fx,Fy,Fz                               ,&
                                   & Fl,Fm,Fn                               ,&
                                   & Fdrag,Fside,Flift                      ,&
                                   & coef_x,coef_y,coef_z                   ,&
                                   & coef_l,coef_m,coef_n                   ,&
                                   & coef_drag,coef_side,coef_lift
! ==============================================================================

contains

! SUBROUTINE PRESSURE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine pressure( interactive                                        ,&
                       & panel_num                                          ,&
                       & panel_num_no_wake                                  ,&
                       & panel_type                                         ,&
                       & case_num                                           ,&
                       & alfa,beta                                          ,&
                       & mac                                                ,&
                       & wing_span                                          ,&
                       & wing_surf                                          ,&
                       & origin                                             ,&
                       & speed                                              ,&
                       & mach                                               ,&
                       & ro                                                 ,&
                       & p_ref                                              ,&
                       & S                                                  ,&
                       & v                                                  ,&
                       & n1,n2,n3                                           ,&
                       & cx,cy,cz                                           )
    
    implicit none
    
    ! INTENTS IN ===================================================================
    integer,                                intent(in) :: panel_num         ,&
                                                        & panel_num_no_wake ,&
                                                        & interactive       ,&
                                                        & case_num
    integer, dimension(panel_num),          intent(in) :: panel_type
    real,                                   intent(in) :: mac               ,&
                                                        & wing_span         ,&
                                                        & wing_surf         ,&
                                                        & speed             ,&
                                                        & mach              ,&
                                                        & ro                ,&
                                                        & p_ref
    real,    dimension(panel_num),          intent(in) :: S                 ,&
                                                        & n1,n2,n3          ,&
                                                        & cx,cy,cz
    real,    dimension(case_num),           intent(in) :: alfa,beta
    real,    dimension(3),                  intent(in) :: origin
    real,    dimension(panel_num_no_wake,&
                      &case_num),           intent(in) :: v
    ! PRIVATE ======================================================================
    integer                                            :: i,j               ,&
                                                        & alloc_stat
    real                                               :: q                 ,&
                                                        & dx,dy,dz          ,&
                                                        & v_square          ,&
                                                        & speed_square      ,&
                                                        & beta_mach
    ! ==============================================================================
    
    speed_square=speed**2
    q=0.5*ro*speed_square
    if (mach .gt. 0.) then
        beta_mach=sqrt(abs(1.-mach**2))
    endif
    
    allocate (     cp(panel_num_no_wake,case_num)   ,&
             & p_dyna(panel_num_no_wake,case_num)   ,&
             & p_mano(panel_num_no_wake,case_num)   ,&
             & p_stat(panel_num_no_wake,case_num)   ,&
             & Fx(case_num)                         ,&
             & Fy(case_num)                         ,&
             & Fz(case_num)                         ,&
             & Fl(case_num)                         ,&
             & Fm(case_num)                         ,&
             & Fn(case_num)                         ,&
             & Fdrag(case_num)                      ,&
             & Fside(case_num)                      ,&
             & Flift(case_num)                      ,&
             & coef_x(case_num)                     ,&
             & coef_y(case_num)                     ,&
             & coef_z(case_num)                     ,&
             & coef_l(case_num)                     ,&
             & coef_m(case_num)                     ,&
             & coef_n(case_num)                     ,&
             & coef_drag(case_num)                  ,&
             & coef_side(case_num)                  ,&
             & coef_lift(case_num)                  ,&
             & stat=alloc_stat                      )
    if (alloc_stat .ne. 0) then
        call func_message( interactive                                          ,&
                         & "    ERROR: Not enough memory for allocating pressure and force arrays")
        press_err=1
        goto 999
    else
        press_err=0
    endif
    
    do j=1,case_num
        Fx(j)=0.
        Fy(j)=0.
        Fz(j)=0.
        Fl(j)=0.
        Fm(j)=0.
        Fn(j)=0.
        do i=1,panel_num_no_wake
            ! avoid dummy panels
            if (panel_type(i) .ne. 20 .and. panel_type(i) .ne. 21) then
                v_square=v(i,j)**2
                
                ! calculating pressure coefficient
                if (mach .gt. 0.) then
                    cp(i,j)=(1.-v_square/speed_square)/beta_mach
                else
                    cp(i,j)=1.-v_square/speed_square
                endif
                
                ! calculating pressures
                p_dyna(i,j)=0.5*ro*v_square
                p_mano(i,j)=cp(i,j)*q
                p_stat(i,j)=p_mano(i,j)+p_ref
                
                ! calculating force components
                dx=-p_mano(i,j)*S(i)*n1(i)
                dy=-p_mano(i,j)*S(i)*n2(i)
                dz=-p_mano(i,j)*S(i)*n3(i)
                
                ! summing force components
                Fx(j)=Fx(j) + dx
                Fy(j)=Fy(j) + dy
                Fz(j)=Fz(j) + dz
                
                Fl(j)=Fl(j)                        - dy*(cz(i)-origin(3)) + dz*(cy(i)-origin(2))
                Fm(j)=Fm(j) + dx*(cz(i)-origin(3))                        - dz*(cx(i)-origin(1))
                Fn(j)=Fn(j) - dx*(cy(i)-origin(2)) + dy*(cx(i)-origin(1))
            endif
        enddo
        
        ! calculating coefficients in body reference frame
        coef_x(j)=Fx(j)/q/wing_surf
        coef_y(j)=Fy(j)/q/wing_surf
        coef_z(j)=Fz(j)/q/wing_surf
        coef_l(j)=Fl(j)/q/wing_surf/wing_span
        coef_m(j)=Fm(j)/q/wing_surf/mac
        coef_n(j)=Fn(j)/q/wing_surf/wing_span
        
        ! calculating forces and coefficients in aerodynamic reference frame
        Fdrag(j)=  Fx(j)*cos(alfa(j))*cos(beta(j)) + Fy(j)*sin(beta(j)) + Fz(j)*sin(alfa(j))*cos(beta(j))
        Fside(j)=  Fx(j)*cos(alfa(j))*sin(beta(j)) + Fy(j)*cos(beta(j)) + Fz(j)*sin(alfa(j))*sin(beta(j))
        Flift(j)=- Fx(j)*sin(alfa(j))                                   + Fz(j)*cos(alfa(j))
        coef_drag(j)=Fdrag(j)/q/wing_surf
        coef_side(j)=Fside(j)/q/wing_surf
        coef_lift(j)=Flift(j)/q/wing_surf
    enddo
        
999 continue
    
    return
    end subroutine pressure

end module module_pressure
