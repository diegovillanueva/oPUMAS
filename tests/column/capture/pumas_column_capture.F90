! 🔬📥 Records the arguments of every process call inside micro_pumas_tend.
! 🧷📄 patch_v1.py inserts capture() before (phase 'in') and after (phase 'out') each call;
!     the file named by $PUMAS_CAPTURE_FILE gets one variable <site>__<arg>__<phase> per argument.
! 🧊📏 Hydrometeor property structs are split into their fields by patch_v1.py.
! 🎛️🔕 Only the capture build of the column driver links this; CESM never sees it.
module pumas_column_capture
  use pumas_kinds, only: r8 => kind_r8
  use netcdf
  implicit none
  private
  public :: capture, capture_site, capture_close

  interface capture
    module procedure cap_real, cap_int, cap_logical
  end interface capture

  integer, save :: ncid = -1, mgncol = 0, nlev = 0
  logical, save :: opened = .false., disabled = .false.
  character(len=64), save :: site = ''
  character(len=64), save :: visited(400) = ''
  integer, save :: nvisits(400) = 0, nsites = 0

contains

  subroutine check(status, what)
    integer, intent(in) :: status
    character(len=*), intent(in) :: what
    if (status /= nf90_noerr) then
      write(*,*) 'capture: ', trim(what), ': ', trim(nf90_strerror(status))
      error stop 1
    end if
  end subroutine check

  ! 🏷️🔢 Start a call site; the second visit of the same site gets suffix _v2, and so on.
  subroutine capture_site(name, ncol, nlev_in, phase)
    character(len=*), intent(in) :: name, phase
    integer, intent(in) :: ncol, nlev_in
    character(len=512) :: path
    integer :: i, k, st
    if (disabled) return
    if (.not. opened) then
      call get_environment_variable('PUMAS_CAPTURE_FILE', path, status=st)
      if (st /= 0 .or. len_trim(path) == 0) then
        disabled = .true.
        return
      end if
      call check(nf90_create(trim(path), ior(nf90_clobber, nf90_netcdf4), ncid), path)
      mgncol = ncol
      nlev = nlev_in
      call check(nf90_put_att(ncid, nf90_global, 'nlev', nlev), 'nlev')
      call check(nf90_enddef(ncid), 'enddef')
      opened = .true.
    end if
    k = 0
    do i = 1, nsites
      if (visited(i) == name) k = i
    end do
    if (k == 0) then
      nsites = nsites + 1
      k = nsites
      visited(k) = name
    end if
    if (phase == 'in') nvisits(k) = nvisits(k) + 1
    if (nvisits(k) == 1) then
      site = name
    else
      write(site, '(a,"_v",i0)') trim(name), nvisits(k)
    end if
  end subroutine capture_site

  subroutine capture_close()
    if (opened) call check(nf90_close(ncid), 'close')
    opened = .false.
  end subroutine capture_close

  ! 📐🗂️ Dimension ids for a shape: col, lev, ilev for model shapes, d<n> otherwise.
  subroutine dims_for(shp, dimids)
    integer, intent(in) :: shp(:)
    integer, allocatable, intent(out) :: dimids(:)
    integer :: i, n
    character(len=32) :: dname
    n = size(shp)
    if (n == 1) then
      ! 🧱↔️ A 1-D vector over all packed points is stored as (col, lev).
      if (shp(1) == mgncol * nlev) then
        allocate(dimids(2))
        call dim_id('col', mgncol, dimids(1))
        call dim_id('lev', nlev, dimids(2))
        return
      end if
    end if
    allocate(dimids(n))
    do i = 1, n
      if (i == 1 .and. shp(i) == mgncol) then
        dname = 'col'
      else if (shp(i) == nlev) then
        dname = 'lev'
      else if (shp(i) == nlev + 1) then
        dname = 'ilev'
      else
        write(dname, '("d",i0)') shp(i)
      end if
      call dim_id(trim(dname), shp(i), dimids(i))
    end do
  end subroutine dims_for

  subroutine dim_id(name, length, id)
    character(len=*), intent(in) :: name
    integer, intent(in) :: length
    integer, intent(out) :: id
    if (nf90_inq_dimid(ncid, name, id) /= nf90_noerr) then
      call check(nf90_redef(ncid), 'redef')
      call check(nf90_def_dim(ncid, name, length, id), name)
      call check(nf90_enddef(ncid), 'enddef')
    end if
  end subroutine dim_id

  integer function new_var(key, phase, xtype, shp, expr) result(varid)
    character(len=*), intent(in) :: key, phase, expr
    integer, intent(in) :: xtype, shp(:)
    integer, allocatable :: dimids(:)
    character(len=200) :: name
    name = trim(site)//'__'//trim(key)//'__'//trim(phase)
    call dims_for(shp, dimids)
    call check(nf90_redef(ncid), 'redef')
    if (size(shp) == 0) then
      call check(nf90_def_var(ncid, trim(name), xtype, varid), name)
    else
      call check(nf90_def_var(ncid, trim(name), xtype, dimids, varid), name)
    end if
    call check(nf90_put_att(ncid, varid, 'expr', trim(expr)), 'expr')
    call check(nf90_enddef(ncid), 'enddef')
  end function new_var

  subroutine cap_real(key, phase, expr, x)
    character(len=*), intent(in) :: key, phase, expr
    real(r8), intent(in) :: x(..)
    integer :: varid
    if (.not. opened) return
    varid = new_var(key, phase, nf90_double, shape(x), expr)
    select rank (x)
    rank (0)
      call check(nf90_put_var(ncid, varid, x), key)
    rank (1)
      if (size(x) == mgncol * nlev) then
        call check(nf90_put_var(ncid, varid, reshape(x, [mgncol, nlev])), key)
      else
        call check(nf90_put_var(ncid, varid, x), key)
      end if
    rank (2)
      call check(nf90_put_var(ncid, varid, x), key)
    rank default
      error stop 'capture: rank > 2'
    end select
  end subroutine cap_real

  subroutine cap_int(key, phase, expr, x)
    character(len=*), intent(in) :: key, phase, expr
    integer, intent(in) :: x(..)
    integer :: varid
    if (.not. opened) return
    varid = new_var(key, phase, nf90_int, shape(x), expr)
    select rank (x)
    rank (0)
      call check(nf90_put_var(ncid, varid, x), key)
    rank (1)
      call check(nf90_put_var(ncid, varid, x), key)
    rank default
      error stop 'capture: integer rank > 1'
    end select
  end subroutine cap_int

  subroutine cap_logical(key, phase, expr, x)
    character(len=*), intent(in) :: key, phase, expr
    logical, intent(in) :: x
    integer :: varid
    if (.not. opened) return
    varid = new_var(key, phase, nf90_int, [integer ::], expr)
    call check(nf90_put_var(ncid, varid, merge(1, 0, x)), key)
  end subroutine cap_logical

end module pumas_column_capture
