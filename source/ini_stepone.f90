subroutine stepone()
    ! subroutine stepone
    !
    ! purpose : call initialization of semi-implicit scheme
    !           and perform initial time step

    use mod_tsteps, only: delt, delt2, alph, rob, wil
    use mod_date, only: istart
    use rp_emulator
    use mod_prec, only: dp

    implicit none

    type(rpe_var) :: delth

    if (istart==0 .or. istart==2) then

      delth = 0.5_dp * delt

      ! semi-impl. initialization
      call impint(delth, alph)

      ! forward half-step
      call step(1, 1, delth, alph, rob, wil)

      ! semi-impl. initialization
      call impint(delt, alph)

      ! leapfrog half-step
      call step(1, 2, delt, alph, rob, wil)
    end if

    ! semi-impl. initialization
    call impint(delt2, alph)
end subroutine stepone
