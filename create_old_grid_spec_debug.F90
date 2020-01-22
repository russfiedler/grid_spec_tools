program create_old_grid_spec
!
! Create a file with all the wet points that we wish to map to.
!
! Mosaic input.
!
! Usage: create_old_grid_spec
!
! Assumes there's a file 'mosaic.nc' lurking about....
!
! Output. A file called old_grid_spec.nc
! At the moment I only calculate the main locations but not the distances
! ds_XX_YY which
! can easily be calculated fronm the folloing table which
! appears in ocean_grids.F90
!
! Also need the angle, areas and Wet.
!
! Note that not all the ds_XX_YY and areas in a full grid_spec.nc are required by MOM however they may be required in generating exchange grids etc.
! It should be straight forward to add these.
!
!---------------------------------------------------------------------------------------------------------------------
!
! Define horizontal (and some vertical) grid arrays.
!
!---------------------------------------------------------------------------------------------------------------------
! Grid%       grid_spec      grid_spec     grid_spec                       Description
!  var        field          field         field
!             VERSION_0      VERSION_1     VERSION_2(mosaic) 
!--------------------------------------------------------------------------------------------------------------------- 
!                                                    ocean_vgrid.nc
!                                                    k=1,nk
! zt          zt             zt            zeta(2k-1)
! zw          zw             zb            zeta(2k)
!
!                                                    ocean_hgrid.nc
!                                                    i=1,ni
!                                                    j=1,nj
! grid_x_t    gridlon_t      grid_x_T      x(2i  ,2)
! grid_x_u    gridlon_vert_t grid_x_C      x(2i+1,1)
! grid_y_t    gridlat_t      grid_y_T      y(ni/4,2j)
! grid_y_u    gridlat_vert_t grid_y_C      y(ni/4,2j+1)
!
!T
! xt(i,j)     geolon_t(i,j)  x_T(i,j)      x(2i,2j)
! yt          geolat_t       y_T           y(2i,2j)
! dtw         dtw            ds_01_11_T    dx(2i-1,2j)                 distance to western face of t cell
! dte         dte            ds_11_21_T    dx(2i,2j)                   distance to eastern face of t cell
! dts         dts            ds_10_11_T    dy(2i,2j-1)                 distance to southern face of t cell 
! dtn         dtn            ds_11_12_T    dy(2i,2j)                   distance to northern face of t cell
! 
! dxt         dxt            ds_01_21_T    dx(2i,2j)    +dx(2i-1,2j)   width of t cell   
! dxtn        dxtn           ds_02_22_T    dx(2i-1,2j+1)+dx(2i,2j+1)   width of northern face of t cell
! dxte        dxte           ds_00_20_C    dx(2i,2j)    +dx(2i+1,2j)   distance to adjacent t cell to the east!
! dyt         dyt            ds_10_12_T    dy(2i,2j)    +dy(2i,2j-1)   height of t cell
! dytn        dytn           ds_00_02_C    dy(2i,2j)    +dy(2i,2j+1)   distance to adjacent t cell to the north!
! dyte        dyte           ds_20_22_T    dy(2i+1,2j-1)+dy(2i+1,2j)   height of eastern face of t cell 
!
!C 
! NOTE: The "first" (I,J) C-cell is the one shifted NE of the "first" (I,J) T-cell
!            
!
! xu          geolon_c       x_C           x(2i+1,2j+1)
! yu          geolat_c       y_c           y(2i+1,2j+1)
! dxu         dxu            ds_01_21_C    dx(2i+1,2j+1)+dx(2i,2j+1)   width of u cell
! dxun        dxun           ds_02_22_C    dx(2i,2j+2)+dx(2i+1,2j+2)   width of northern face of u cell
! dyu         dyu            ds_10_12_C    dy(2i+1,2j+1)+dy(2i+1,2j)   height of u cell
! dyue        dyue           ds_20_22_C    dy(2i+2,2j)+dy(2i+2,2j+1)   height of eastern face of u cell
!
! dyun        dyun           ds_11_12_C    dy(2i+1,2j+1)+dy(2i+1,2j+2) distance to adjacent u cell to the north  !RASF see below 
!                            +ds_10_11_C(i,j+1)                         satisfies sum rule dyte(i,j)=dyun(i,j-1) 
! dxue        dxue           ds_11_21_C    dx(2i+1,2j+1)+dx(2i+2,2j+1) distance to adjacent u cell to the east!  ! RASF see below
!                            +ds_01_11_C(i+1,j)   
!
! duw         duw            ds_01_11_C    dx(2i,2j+1)                 distance to western face of u cell
! due         due            ds_11_21_C    dx(2i+1,2j+1)               distance to eastern face of u cell
! dus         dus            ds_10_11_C    dy(2i+1,2j)                 distance to southern face of u cell
! dun         dun            ds_11_12_C    dy(2i+1,2j+1)               distance to northern face of u cell
!  
! sin_rot     sin_rot        angle_C       sin(angle_dx(2*i+1,2*j+1)   sin of rotation angle at corner cell centers
! cos_rot     cos_rot        angle_C       cos(angle_dx(2*i+1,2*j+1)   cos of rotation angle at corner cell centers
!
!Following are the available fields in mosaic files
!--------------------------------------------------------
!Mosaic file     fields
!--------------------------------------------------------
!ocean_hgrid.nc  x, y, dx, dy, angle_dx, area
!ocean_vgrid.nc  zeta
!topog.nc        depth
!
!
! Halo updates. Some of the distance quantities require halo information since
! dx is dimensioned (nx,ny+1) and dy is dimensioned (nx+1,ny). i.e. we need to
! add a right hand condition for dx and a top condition for dy.
!
! For tripolar grids we usually have (by symmetry) dy(:,ny+1)=dy(:,ny),
! (imore strictly dy(:,ny+1)=dy(nx+1:1,:-1,ny))
! Cylic condition dx(nx+1,:)=dx(1,:)
! 
!-----------------------------------------------------------------------------------------------------------
!
!subroutine set_ocean_hgrid_arrays(Domain, Grid)

!
!
use iso_fortran_env
use netcdf
implicit none

type grid_type
   integer                               :: nx,ny,nz
   real(real64),allocatable,dimension(:) :: zt
   real(real64),allocatable,dimension(:) :: zb
   real(real64),allocatable,dimension(:,:) :: xt
   real(real64),allocatable,dimension(:,:) :: xc
   real(real64),allocatable,dimension(:,:) :: yt
   real(real64),allocatable,dimension(:,:) :: yc
   real(real64),allocatable,dimension(:,:) :: depth_t
   real(real64),allocatable,dimension(:,:) :: num_levels
end type grid_type

type(grid_type)     :: grid

integer(kind=int32) :: i,j,k
integer(kind=int32) :: nx,ny,nz                  ! Size of model grid
integer(kind=int32) :: nxp,nyp,nzp                ! Size of model supergrid

integer(kind=int32) :: ncid,vid,did           ! NetCDF ids

real(kind=real64),allocatable,dimension(:,:)   :: wrk_super
real(kind=real64),allocatable,dimension(:)   :: zeta

character(len=128)  :: odir='',ofile=''  ! for ocean_mosaic.nc
character(len=128)  :: gdir='',gfile=''  ! for hgrid file
character(len=128)  :: tdir='',tfile=''  ! for topgraphy file
character(len=256)  :: dirfile=''        ! concatenation

logical             :: fexist = .false.

real(real64) :: tol=0.0

! read a tolerance 
write(*,*) 'Enter tolerance'
read(*,*) tol

! Get info on the grid from input

write(*,*) 'Getting model grid info'
! Get mosaic info
inquire(file=trim('mosaic.nc'),exist=fexist)
if ( .not. fexist ) then
   write(*,*) 'mosaic.nc does not exist. Bailing out' 
   stop 1
endif
call handle_error(nf90_open('mosaic.nc',nf90_nowrite,ncid))
call handle_error(nf90_inq_varid(ncid,'ocn_mosaic_dir',vid))
call handle_error(nf90_get_var(ncid,vid,odir))
call handle_error(nf90_inq_varid(ncid,'ocn_mosaic_file',vid))
call handle_error(nf90_get_var(ncid,vid,ofile))
call handle_error(nf90_inq_varid(ncid,'ocn_topog_dir',vid))
call handle_error(nf90_get_var(ncid,vid,tdir))
call handle_error(nf90_inq_varid(ncid,'ocn_topog_file',vid))
call handle_error(nf90_get_var(ncid,vid,tfile))
call handle_error(nf90_close(ncid))
! Get horizontal grid
dirfile=odir(1:scan(odir,'/',back=.true.)) // ofile(1:scan(ofile,'c',back=.true.))
write(*,*) len_trim(dirfile),dirfile
inquire(file=trim(dirfile),exist=fexist)
if ( .not. fexist ) then
   write(*,*) 'ocn_mosaic_dir/ocn_mosaic_file =',trim(dirfile), ' does not exist. Bailing out' 
   stop 1
endif
call handle_error(nf90_open(trim(dirfile),nf90_nowrite,ncid))
call handle_error(nf90_inq_varid(ncid,'gridlocation',vid))
call handle_error(nf90_get_var(ncid,vid,gdir))
call handle_error(nf90_inq_varid(ncid,'gridfiles',vid))
call handle_error(nf90_get_var(ncid,vid,gfile))
call handle_error(nf90_close(ncid))


dirfile=tdir(1:scan(tdir,'/',back=.true.)) // tfile(1:scan(tfile,'c',back=.true.))
inquire(file=trim(dirfile),exist=fexist)
if ( .not. fexist ) then
   write(*,*) 'ocn_topog_dir/ocn_topog_file =',trim(dirfile), ' does not exist. Bailing out' 
   stop 1
endif
call handle_error(nf90_open(trim(dirfile),nf90_nowrite,ncid))
! Kial has a grid file with 'xx' and 'yy' for dimensions and variables
call handle_error(nf90_inq_dimid(ncid,'xx',did))
call handle_error(nf90_inquire_dimension(ncid,did,len=nx))
call handle_error(nf90_inq_dimid(ncid,'yy',did))
call handle_error(nf90_inquire_dimension(ncid,did,len=ny))
call handle_error(nf90_inq_varid(ncid,'depth',vid))

allocate(grid%xt(nx,ny),grid%xc(nx,ny),grid%yt(nx,ny),grid%yc(nx,ny),grid%depth_t(nx,ny),grid%num_levels(nx,ny))

! Get depth

call handle_error(nf90_get_var(ncid,vid,grid%depth_t))
call handle_error(nf90_close(ncid))

!
! On mosaic "supergrid" we need to get every second point
!
write(*,*) 'Reading supergrid info'
nxp = 2*nx+1
nyp = 2*ny+1
allocate(wrk_super(nxp,nyp))
! Read xt
dirfile=gdir(1:scan(gdir,'/',back=.true.)) // gfile(1:scan(gfile,'c',back=.true.))
inquire(file=trim(dirfile),exist=fexist)
if ( .not. fexist ) then
   write(*,*) 'gridlocation/gridfiles =',trim(dirfile), ' does not exist. Bailing out' 
   stop 1
endif


! Get x points

call handle_error(nf90_open(trim(dirfile),nf90_nowrite,ncid))
call handle_error(nf90_inq_varid(ncid,'x',vid))
call handle_error(nf90_get_var(ncid,vid,wrk_super))
do j=1,ny
   do i = 1,nx
      grid%xt(i,j)= wrk_super(2*i,2*j)
   enddo
enddo
do j=1,ny
   do i = 1,nx
      grid%xc(i,j)= wrk_super(2*i+1,2*j+1)
   enddo
enddo

! Handle y values
!
call handle_error(nf90_inq_varid(ncid,'y',vid))
call handle_error(nf90_get_var(ncid,vid,wrk_super))
do j=1,ny
   do i = 1,nx
      grid%yt(i,j)= wrk_super(2*i,2*j)
   enddo
enddo
do j=1,ny
   do i = 1,nx
      grid%yc(i,j)= wrk_super(2*i+1,2*j+1)
   enddo
enddo
call handle_error(nf90_close(ncid))



!
! Get zt
!

! Read xt
dirfile=gdir(1:scan(gdir,'/',back=.true.)) // 'ocean_vgrid.nc'
inquire(file=trim(dirfile),exist=fexist)
if ( .not. fexist ) then
   write(*,*) 'gridlocation/gridfiles =',trim(dirfile), ' does not exist. Bailing out' 
   stop 1
endif
call handle_error(nf90_open(trim(dirfile),nf90_nowrite,ncid))
call handle_error(nf90_inq_dimid(ncid,'nzv',did))
call handle_error(nf90_inquire_dimension(ncid,did,len=nzp))
nz=(nzp-1)/2
allocate(zeta(nzp),grid%zt(nz),grid%zb(nz))
call handle_error(nf90_inq_varid(ncid,'zeta',vid))
call handle_error(nf90_get_var(ncid,vid,zeta))
call handle_error(nf90_close(ncid))
do k = 1, nz
   grid%zt(k)=zeta(2*k)
   grid%zb(k)=zeta(2*k+1)
enddo

! Now get number of nevels
grid%num_levels=0
do j=1,ny
   do i=1,nx
      if(grid%depth_t(i,j) <= -1.0e6_real64)  grid%depth_t(i,j)=0.0
      if(grid%depth_t(i,j) <= 0.0_real64) cycle
      do k=1,nz
         if(grid%depth_t(i,j) <= grid%zb(k)+tol) then
            grid%num_levels(i,j)=k
            exit
         endif
      enddo
   enddo
enddo
            
grid%nx=nx
grid%ny=ny
grid%nz=nz

write(*,*) 'Writing'
call create_grid_spec_file(grid)
write(*,*) 'Done'

contains

subroutine create_grid_spec_file(grid)
   type(grid_type), intent(in) :: grid
   integer(kind=int32) :: ncid
   integer(kind=int32) :: did_xt,did_xc,did_yt,did_yc,did_zt,did_zb
   integer(kind=int32) :: zt_id,zb_id
   integer(kind=int32) :: xt_id,xc_id
   integer(kind=int32) :: yt_id,yc_id
   integer(kind=int32) :: nl_id,dt_id


!   call handle_error(nf90_create('old_grid_spec.nc',ior(NF90_CLOBBER,NF90_HDF5),ncid))  ! Deprecated. For old versions of netCDF only.
   call handle_error(nf90_create('old_grid_spec.nc',ior(NF90_CLOBBER,NF90_NETCDF4),ncid))
   call handle_error(nf90_def_dim(ncid,'grid_x_T',grid%nx,did_xt))
   call handle_error(nf90_def_dim(ncid,'grid_y_T',grid%ny,did_yt))
   call handle_error(nf90_def_dim(ncid,'grid_x_C',grid%nx,did_xc))
   call handle_error(nf90_def_dim(ncid,'grid_y_C',grid%ny,did_yc))
   call handle_error(nf90_def_dim(ncid,'zt',grid%nz,did_zt))
   call handle_error(nf90_def_dim(ncid,'zb',grid%nz,did_zb))
!
!
   call handle_error(nf90_def_var(ncid,'zt',nf90_real,did_zt,zt_id))
   call handle_error(nf90_def_var(ncid,'zb',nf90_real,did_zb,zb_id))

   call handle_error(nf90_def_var(ncid,'x_T',nf90_double,(/did_xt,did_yt/),xt_id))
   call handle_error(nf90_put_att(ncid,xt_id,'long_name','Geographic longitude of T_cell centers'))
   call handle_error(nf90_put_att(ncid,xt_id,'units','degrees_E'))
   call handle_error(nf90_def_var(ncid,'x_C',nf90_double,(/did_xc,did_yc/),xc_id))
   call handle_error(nf90_put_att(ncid,xc_id,'long_name','Geographic longitude of C_cell centers'))
   call handle_error(nf90_put_att(ncid,xc_id,'units','degrees_E'))
   call handle_error(nf90_def_var(ncid,'y_T',nf90_double,(/did_xt,did_yt/),yt_id))
   call handle_error(nf90_put_att(ncid,yt_id,'long_name','Geographic latitude of T_cell centers'))
   call handle_error(nf90_put_att(ncid,yt_id,'units','degrees_N'))
   call handle_error(nf90_def_var(ncid,'y_C',nf90_double,(/did_xc,did_yc/),yc_id))
   call handle_error(nf90_put_att(ncid,yc_id,'long_name','Geographic latitude of C_cell centers'))
   call handle_error(nf90_put_att(ncid,yc_id,'units','degrees_N'))
   call handle_error(nf90_def_var(ncid,'depth_t',nf90_double,(/did_xt,did_yt/),dt_id))
   call handle_error(nf90_put_att(ncid,dt_id,'long_name','topographic depth of T-cell'))
   call handle_error(nf90_put_att(ncid,dt_id,'units','meters'))
   call handle_error(nf90_def_var(ncid,'num_levels',nf90_double,(/did_xt,did_yt/),nl_id))
   call handle_error(nf90_put_att(ncid,nl_id,'long_name','number of vertical T-cells'))
   call handle_error(nf90_put_att(ncid,nl_id,'units','none'))




   call handle_error(nf90_enddef(ncid,h_minfree=4096))

! Put it there
   call handle_error(nf90_put_var(ncid,xt_id,grid%xt))
   call handle_error(nf90_put_var(ncid,xc_id,grid%xc))
   call handle_error(nf90_put_var(ncid,yt_id,grid%yt))
   call handle_error(nf90_put_var(ncid,yc_id,grid%yc))
   call handle_error(nf90_put_var(ncid,nl_id,grid%num_levels))
   call handle_error(nf90_put_var(ncid,dt_id,grid%depth_t))
   call handle_error(nf90_put_var(ncid,zt_id,grid%zt))
   call handle_error(nf90_put_var(ncid,zb_id,grid%zb))

   call handle_error(nf90_close(ncid))

end subroutine create_grid_spec_file

subroutine handle_error(error_flag,isfatal,err_string)
#ifdef __INTEL_COMPILER
!  Needs special traceback module for intel compile RASF
   use ifcore
#endif
! Simple error handle for netCDF
integer(kind=int32),intent(in) :: error_flag
logical, intent(in),optional :: isfatal
character(*), intent(in),optional :: err_string
logical            :: fatal
fatal = .true.
if(present(isfatal)) fatal=isfatal
if ( error_flag  /= nf90_noerr ) then
   if ( fatal ) then
      write(*,*) 'FATAL ERROR:',nf90_strerror(error_flag)
      if (present(err_string)) write(*,*) trim(err_string)
#ifdef __INTEL_COMPILER
        ! Get traceback and return quietly for correct abort
        call TRACEBACKQQ(user_exit_code=-1)
#endif
      stop
   endif
endif
end subroutine handle_error
    
end program create_old_grid_spec
