program main
    use driver
    implicit none
    
    ! so far the only thing I have main doing is passing file names into
    ! the driver module
    character(len=256) :: config_file = 'triso_config.txt'
    character(len=256) :: xs_file_list = 'csv_xs_list.txt' ! now .dat, not .csv
end program main

