!-----------------------------------------------------------------------------
! (c) Crown copyright 2025 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Computes the RHS for projecting orography into W0 function space
!> @details This kernel assembles the right-hand side vector for the projection
!!          of a surface altitude field into W3 into W0, such that the W3
!!          surface altitude field is preserved.
!!          This involves integrating W0 test functions against the W3 field
!!          on the top/bottom faces of a 2D mesh.
!!          Note: this assumes the input surface altitude field is in the
!!          lowest-order W3 space
module w0_orography_combine_kernel_mod

  use argument_mod,             only : arg_type, func_type,                    &
                                       GH_FIELD, GH_REAL, GH_READ, GH_WRITE,   &
                                       GH_BASIS, GH_DIFF_BASIS,                &
                                       GH_EVALUATOR, CELL_COLUMN, ANY_SPACE_9
  use constants_mod,            only : r_def, i_def, l_def, EPS
  use fs_continuity_mod,        only : W0, W3, W1, W2
  use kernel_mod,               only : kernel_type
  use reference_element_mod,    only : W, E, S, N, B, T, WB, EB, SB, NB,       &
                                       SW, SE, NE, NW, SWB, SEB, NEB, NWB,     &
                                       SWT, SET, NET, NWT, ST, NT, WT, ET

  implicit none

  private

  !-------------------------------------------------------------------------------
  ! Public types
  !-------------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the PSy layer
  type, public, extends(kernel_type) :: w0_orography_combine_kernel_type
    private
    type(arg_type) :: meta_args(6) = (/                                        &
        arg_type(GH_FIELD,   GH_REAL,    GH_WRITE, ANY_SPACE_9),               &
        arg_type(GH_FIELD,   GH_REAL,    GH_READ,  W0),                        &
        arg_type(GH_FIELD,   GH_REAL,    GH_READ,  W1),                        &
        arg_type(GH_FIELD,   GH_REAL,    GH_READ,  W2),                        &
        arg_type(GH_FIELD,   GH_REAL,    GH_READ,  W3),                        &
        arg_type(GH_FIELD,   GH_REAL,    GH_READ,  ANY_SPACE_9)                &
    /)
    type(func_type) :: meta_funcs(5) = (/                                      &
        func_type(ANY_SPACE_9, GH_BASIS),                                      &
        func_type(W0, GH_BASIS),                                               &
        func_type(W1, GH_BASIS),                                               &
        func_type(W2, GH_BASIS),                                               &
        func_type(W3, GH_BASIS)                                                &
    /)
    integer :: operates_on = CELL_COLUMN
    integer :: gh_shape = GH_EVALUATOR
  contains
    procedure, nopass :: w0_orography_combine_code
  end type

  !-------------------------------------------------------------------------------
  ! Contained functions/subroutines
  !-------------------------------------------------------------------------------
  public :: w0_orography_combine_code

contains

  !> @brief Computes the RHS for projecting orography into W0 function space
  subroutine w0_orography_combine_code(nlayers,                                &
                                       surface_altitude_w0,                    &
                                       surface_altitude_w0_k0,                 &
                                       surface_altitude_w1_k0,                 &
                                       dummy_w2_k0,                            &
                                       surface_altitude_w3_k0,                 &
                                       surface_altitude_w0_avg,                &
                                       ndf_w0, undf_w0, map_w0,                &
                                       basis_w0,                               &
                                       ndf_w0_k0, undf_w0_k0, map_w0_k0,       &
                                       basis_w0_k0,                            &
                                       ndf_w1_k0, undf_w1_k0, map_w1_k0,       &
                                       basis_w1_k0,                            &
                                       ndf_w2_k0, undf_w2_k0, map_w2_k0,       &
                                       basis_w2_k0,                            &
                                       ndf_w3_k0, undf_w3_k0, map_w3_k0,       &
                                       basis_w3_k0)

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in) :: nlayers
    integer(kind=i_def), intent(in) :: ndf_w0, ndf_w3_k0, ndf_w0_k0, ndf_w1_k0, ndf_w2_k0
    integer(kind=i_def), intent(in) :: undf_w0, undf_w3_k0, undf_w0_k0, undf_w1_k0, undf_w2_k0

    integer(kind=i_def), intent(in) :: map_w0(ndf_w0)
    integer(kind=i_def), intent(in) :: map_w3_k0(ndf_w3_k0)
    integer(kind=i_def), intent(in) :: map_w0_k0(ndf_w0_k0)
    integer(kind=i_def), intent(in) :: map_w1_k0(ndf_w1_k0)
    integer(kind=i_def), intent(in) :: map_w2_k0(ndf_w2_k0)

    real(kind=r_def),    intent(in) :: basis_w0(1,ndf_w0,ndf_w0)
    real(kind=r_def),    intent(in) :: basis_w0_k0(1,ndf_w0_k0,ndf_w0)
    real(kind=r_def),    intent(in) :: basis_w1_k0(3,ndf_w1_k0,ndf_w0)
    real(kind=r_def),    intent(in) :: basis_w2_k0(3,ndf_w2_k0,ndf_w0)
    real(kind=r_def),    intent(in) :: basis_w3_k0(1,ndf_w3_k0,ndf_w0)

    real(kind=r_def),    intent(inout) :: surface_altitude_w0(undf_w0)
    real(kind=r_def),    intent(in)    :: surface_altitude_w0_k0(undf_w0_k0)
    real(kind=r_def),    intent(in)    :: surface_altitude_w1_k0(undf_w1_k0)
    real(kind=r_def),    intent(in)    :: dummy_w2_k0(undf_w2_k0)
    real(kind=r_def),    intent(in)    :: surface_altitude_w3_k0(undf_w3_k0)
    real(kind=r_def),    intent(in)    :: surface_altitude_w0_avg(undf_w0)

    ! Internal variables
    integer(kind=i_def) :: i
    real(kind=r_def)    :: basis_w1_scalar(1,ndf_w1_k0,ndf_w0)
    real(kind=r_def)    :: basis_w2_scalar(1,ndf_w2_k0,ndf_w0)
    logical(kind=l_def) :: found_match
    integer(kind=i_def) :: df_w0, df_w0_k0, df_w1_k0, df_w2_k0

    ! Make scalar W1 basis functions
    basis_w1_scalar(:,:,:) = 0.0_r_def
    basis_w2_scalar(:,:,:) = 0.0_r_def
    do i = 1, 3
      basis_w1_scalar(1,:,:) = basis_w1_scalar(1,:,:) + ABS(basis_w1_k0(i,:,:))
      basis_w2_scalar(1,:,:) = basis_w2_scalar(1,:,:) + ABS(basis_w2_k0(i,:,:))
    end do

    do df_w0 = 1, ndf_w0
      found_match = .false.

      ! We don't know the order of the DoFs, so use the basis functions to
      ! work out which of the input lowest-order W0, W1, W3 DoFs correspond to
      ! the order=1 W0 DoFs

      ! Order 0 W0 DoFs --------------------------------------------------------
      do df_w0_k0 = 1, ndf_w0_k0
        if (ABS(basis_w0_k0(1,df_w0_k0,df_w0) - 1.0_r_def) < EPS) then
          ! To ensure vertical consistency, take only a single W0 for
          ! each vertical edge. For values at top of cell, take the values
          ! from the bottom of the cell
          select case (df_w0_k0)
          case (SWT)
            surface_altitude_w0(map_w0(df_w0)) = surface_altitude_w0_k0(map_w0_k0(SWB))
          case (SET)
            surface_altitude_w0(map_w0(df_w0)) = surface_altitude_w0_k0(map_w0_k0(SEB))
          case (NET)
            surface_altitude_w0(map_w0(df_w0)) = surface_altitude_w0_k0(map_w0_k0(NEB))
          case (NWT)
            surface_altitude_w0(map_w0(df_w0)) = surface_altitude_w0_k0(map_w0_k0(NWB))
          case (SWB, SEB, NEB, NWB)
            ! DoF is already at the bottom, so take this value
            surface_altitude_w0(map_w0(df_w0)) = surface_altitude_w0_k0(map_w0_k0(df_w0_k0))
          end select
          found_match = .true.
          EXIT
        end if
      end do

      if (found_match) CYCLE   ! Skip to the next order=1 W0 DoF

      ! Order 0 W1 DoFs --------------------------------------------------------
      do df_w1_k0 = 1, ndf_w1_k0
        if (ABS(basis_w1_scalar(1,df_w1_k0,df_w0) - 1.0_r_def) < EPS) then
          ! To ensure vertical consistency, ensure that W1 values on vertical
          ! edges match the W0 values. For values horizontal edges, we take
          ! the values from the bottom of the cell
          select case (df_w1_k0)
          ! Vertical edges
          case (SW)
            surface_altitude_w0(map_w0(df_w0)) = surface_altitude_w0_k0(map_w0_k0(SWB))
          case (SE)
            surface_altitude_w0(map_w0(df_w0)) = surface_altitude_w0_k0(map_w0_k0(SEB))
          case (NE)
            surface_altitude_w0(map_w0(df_w0)) = surface_altitude_w0_k0(map_w0_k0(NEB))
          case (NW)
            surface_altitude_w0(map_w0(df_w0)) = surface_altitude_w0_k0(map_w0_k0(NWB))
          ! Horizontal edges
          case (WT)
            surface_altitude_w0(map_w0(df_w0)) = surface_altitude_w1_k0(map_w1_k0(WB))
          case (ET)
            surface_altitude_w0(map_w0(df_w0)) = surface_altitude_w1_k0(map_w1_k0(EB))
          case (ST)
            surface_altitude_w0(map_w0(df_w0)) = surface_altitude_w1_k0(map_w1_k0(SB))
          case (NT)
            surface_altitude_w0(map_w0(df_w0)) = surface_altitude_w1_k0(map_w1_k0(NB))
          case (WB, EB, SB, NB)
            ! DoF is already at the bottom, so take this value
            surface_altitude_w0(map_w0(df_w0)) = surface_altitude_w1_k0(map_w1_k0(df_w1_k0))
          end select
          found_match = .true.
          EXIT
        end if
      end do

      if (found_match) CYCLE   ! Skip to the next order=1 W0 DoF

      ! Order 0 W2 DoFs --------------------------------------------------------
      do df_w2_k0 = 1, ndf_w2_k0
        if (ABS(basis_w2_scalar(1,df_w2_k0,df_w0) - 1.0_r_def) < EPS) then
          found_match = .true.

          ! Now take value from (lowest) W1 or W3 fields
          select case (df_w2_k0)
          case (W)
            surface_altitude_w0(map_w0(df_w0)) = surface_altitude_w1_k0(map_w1_k0(WB))
          case (E)
            surface_altitude_w0(map_w0(df_w0)) = surface_altitude_w1_k0(map_w1_k0(EB))
          case (S)
            surface_altitude_w0(map_w0(df_w0)) = surface_altitude_w1_k0(map_w1_k0(SB))
          case (N)
            surface_altitude_w0(map_w0(df_w0)) = surface_altitude_w1_k0(map_w1_k0(NB))
          case (B, T)
            surface_altitude_w0(map_w0(df_w0)) = surface_altitude_w3_k0(map_w3_k0(1))
          end select
        end if
      end do

      ! Order 0 W3 DoFs --------------------------------------------------------
      if (.not. found_match) then
        ! Remaining DoF must be in the centre, so correspond to order 0 W3
        surface_altitude_w0(map_w0(df_w0)) =                                   &
            surface_altitude_w3_k0(map_w3_k0(1))
      end if
    end do

    do df_w0 = 1, ndf_w0
      if (ABS(surface_altitude_w0_avg(map_w0(df_w0))) < EPS) then
        surface_altitude_w0(map_w0(df_w0)) = 0.0_r_def
      end if
    end do

  end subroutine w0_orography_combine_code

end module w0_orography_combine_kernel_mod
