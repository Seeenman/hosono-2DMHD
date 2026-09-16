module assert
    ! module:       assert
    ! purpose:      Minimal, dependency-free assertion helpers for the hand-rolled
    !               unit test programs under test/. Each test program is a small
    !               standalone Fortran `program` that calls into these routines
    !               and finishes by calling assert_summary to report a pass/fail
    !               count and set a process exit status (nonzero on failure).

    implicit none

    private

    public :: assert_true
    public :: assert_equal_int
    public :: assert_close
    public :: assert_summary

    integer :: nChecks = 0
    integer :: nFailed = 0

contains

    subroutine assert_true(condition, message)
        implicit none
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message

        nChecks = nChecks + 1
        if (condition) then
            write(*,'(A,A)') "  [PASS] ", trim(message)
        else
            nFailed = nFailed + 1
            write(*,'(A,A)') "  [FAIL] ", trim(message)
        end if
    end subroutine assert_true

    subroutine assert_equal_int(actual, expected, message)
        implicit none
        integer, intent(in) :: actual, expected
        character(len=*), intent(in) :: message
        character(len=300) :: fullMessage

        write(fullMessage,'(A,A,I0,A,I0,A)') trim(message), " (expected ", expected, ", got ", actual, ")"
        call assert_true(actual==expected, trim(fullMessage))
    end subroutine assert_equal_int

    subroutine assert_close(actual, expected, tol, message)
        implicit none
        real, intent(in) :: actual, expected, tol
        character(len=*), intent(in) :: message
        character(len=300) :: fullMessage

        write(fullMessage,'(A,A,ES13.5,A,ES13.5,A)') &
            trim(message), " (expected ", expected, ", got ", actual, ")"
        call assert_true(ABS(actual-expected)<=tol, trim(fullMessage))
    end subroutine assert_close

    subroutine assert_summary(exitStat)
        implicit none
        integer, intent(out) :: exitStat

        write(*,*) "--------------------------------------------------------------"
        write(*,'(A,I0,A,I0,A)') "Passed ", nChecks-nFailed, " / ", nChecks, " checks"
        if (nFailed>0) then
            write(*,*) "RESULT: FAIL"
            exitStat = 1
        else
            write(*,*) "RESULT: PASS"
            exitStat = 0
        end if
    end subroutine assert_summary

end module assert
