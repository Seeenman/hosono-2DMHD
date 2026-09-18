program test_template
    ! program:      test_template

    use definitions, only: max_string_length
    use assert, only: assert_close, assert_summary

    implicit none

    character(len=max_string_length) :: paramfile

    paramfile = "par/test_template.par"

    call assert_summary(exitStat)
    call EXIT(exitStat)

end program test_template
