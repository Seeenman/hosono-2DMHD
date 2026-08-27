module readParamFile
    ! Module with subroutines that are able to read in different data types from a 
    ! text file. Examples of acceptable ways for each line to 
    ! be formatted are:
    !   # comments
    !   var_name1 =3.2 # a real variable
    !   var_name2= 3 # an integer variable
    !   var_name3 = 'a_string' # a string variable
    !   var_name4 = .true. # a logical (boolean) variable

    use definitions, only: max_string_length

    implicit none

    private

    public :: readParamFile_int 
    public :: readParamFile_real 
    public :: readParamFile_char 
    public :: readParamFile_logical

contains

    function readParamFile_int(file_name, var_name) result(var_value)
        ! function:     readParamFile_int
        ! Author:       Sean Riedel
        ! purpose:      To read in a variable from a .init file
        ! 
        ! Inputs:       - file_name (character) name of file to read
        !               - var_name (character) name of variable as written in .init file
        !               
        ! Outputs:      - var_value (integer) value of variable with name var_name
        ! ------------------------------------------------------------
        implicit none
        character(len=*), intent(in) :: file_name, var_name
        integer :: var_value
        ! local variables
        character(len=max_string_length) :: var_value_char
        
        call read_value_string(file_name, var_name, var_value_char)
        read(var_value_char,*) var_value
    
    end function readParamFile_int

    function readParamFile_real(file_name, var_name) result(var_value)
        ! function:     readParamFile_real
        ! Author:       Sean Riedel
        ! purpose:      To read in a real type variable from a text file
        ! 
        ! Inputs:       - file_name (character) name of file to read
        !               - var_name (character) name of variable as written in .init file
        !               
        ! Outputs:      - var_value (double) value of variable with name var_name
        ! ------------------------------------------------------------
        implicit none
        character(len=*), intent(in) :: file_name, var_name
        real :: var_value
        ! local variables
        character(len=max_string_length) :: var_value_char
        
        call read_value_string(file_name, var_name, var_value_char)
        read(var_value_char,*) var_value

    end function readParamFile_real
    
    function readParamFile_char(file_name, var_name) result(var_value)
        ! function:   readParamFile_char
        ! Author:       Sean Riedel
        ! purpose:      To read in a variable from a .init file
        ! 
        ! Inputs:       - file_name (character) name of file to read
        !               - var_name (character) name of variable as written in .init file
        !               
        ! Outputs:      - var_value (character) value of variable with name var_name
        ! ------------------------------------------------------------
        implicit none
        character(len=*), intent(in) :: file_name, var_name
        character(len=max_string_length) :: var_value
        ! local variables
        character(len=max_string_length) :: var_value_char
        
        call read_value_string(file_name, var_name, var_value_char)
        read(var_value_char,*) var_value
    
    end function readParamFile_char

    function readParamFile_logical(file_name, var_name) result(var_value)
        ! function:   readParamFile_logical
        ! Author:       Sean Riedel
        ! purpose:      To read in a variable from a .init file
        ! 
        ! Inputs:       - file_name (character) name of file to read
        !               - var_name (character) name of variable as written in .init file
        !               
        ! Outputs:      - var_value (logical) value of variable with name var_name
        ! ------------------------------------------------------------
        implicit none
        character(len=*), intent(in) :: file_name, var_name
        logical :: var_value
        ! local variables
        character(len=max_string_length) :: var_value_char
        
        call read_value_string(file_name, var_name, var_value_char)
        read(var_value_char,*) var_value
    
    end function readParamFile_logical

    subroutine read_value_string(file_name, var_name, var_value_char)
        ! subroutine:   read_value_string
        ! Author:       Sean Riedel
        ! purpose:      to return the string that contains the relevant variable from 
        !               an init file
        ! 
        ! Inputs:       - file_name (character) name of file to read
        !               - var_name (character) name of variable as written in .init file
        !               
        ! Outputs:      - var_value_char (character) string containing the value of the
        !                 variable specified by var_name
        ! ------------------------------------------------------------
        implicit none
        character(len=*), intent(in) :: var_name, file_name
        character(len=*), intent(out) :: var_value_char
        ! local variables
        integer :: input_status, open_status
        character(len=max_string_length) :: dummy_var_name, line

        open(unit=10, file=file_name, status='old', IOSTAT=open_status, FORM='formatted', ACTION='read')
        do
            read(10, FMT = 100, IOSTAT = input_status) line ! read the current line of the file into the variable line
            if (line(1:1)/='#') then ! if first character is not a comment character
                if (index(line,'#') > 0) then ! if there is another comment character somewhere else
                    line = trim(line(1:index(line,"#")-1)) ! redefine line to exclude the portion after the #
                end if
                dummy_var_name = trim(line(1:index(line,'=')-1))
                if (dummy_var_name==var_name) then 
                    var_value_char = line(index(line,'=')+1:)
                    exit
                end if
            end if
            if (input_status<0) then ! if the parameter wasn't read in
                write(*,*) "====================================================================================="
                write(*,*) "UNSUCCESFUL INPUT FILE READ." 
                write(*,*) "The parameter            '", trim(var_name), "'"
                write(*,*) "was not found inside of  '", trim(file_name), "'"
                write(*,*) "---------------------------------------------------------------------------------"
                write(*,*) "Make sure the parameter is present in the input"
                write(*,*) "file and is spelled correctly."
                write(*,*) "====================================================================================="
                write(*,*) " "
                stop
            end if
        end do

        write(*,*) trim(var_name)//" set to "//trim(var_value_char)

        close(10)
        100 format(A)
        
    end subroutine read_value_string
        
end module readParamFile
