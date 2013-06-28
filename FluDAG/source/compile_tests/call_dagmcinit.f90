program main
implicit none
character(len=8) :: filename
integer :: clen
character(len=8) :: ftol
integer :: ftlen
integer :: parallel_read
double precision :: dagmc_version
integer :: moab_version
integer :: max_pbl


call dagmcinit(filename,clen,ftol,ftlen,parallel_read,dagmc_version,moab_version,max_pbl)

stop
end program
