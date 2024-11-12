module QuadGrid
  use Globals
  use AbstractGeometry
  implicit none
  private
  !> structure variable containing mesh and geometrical information
  type,  public :: quadmesh
     type(abs_simplex_mesh) :: grid
     !> Number of x interval. **Input in the constructor**.
     integer :: ninter = 0
     !> lower left corner of the square
     real(kind=double) :: base_point(2)
     !> measures x,y of the square
     real(kind=double) :: measures(2) 
     !> Dimension(ninter+1)
     !> x coord. of the nodes
     real(kind=double), allocatable :: xcoord(:)
     !> Dimension(ninter+1)
     !> y coord of the nodes
     real(kind=double), allocatable :: ycoord(:)    
     !> **Input in the constructor.**
     integer, allocatable :: topol_inter(:,:)
   contains
     !> static constructor
     !> (procedure public for type quadmesh)
     procedure, public, pass :: init => init_quadmesh
     !> static destructor
     !> (procedure public for type quadmesh)
     procedure, public, pass :: kill => kill_quadmesh
     !> Procedure for getting x y intervals of a cell
     !> (public procedure for type quadmesh)
     procedure, public, pass :: get_interval_from_cell
  end type quadmesh
contains
  !>-----------------------------------------------
  !> Static destructor.
  !> (procedure public for type quadmesh)
  !> deallocate all arrays for a var of type quadmesh
  !>
  !> usage:
  !>     call 'var'%kill(lun)
  !>
  !> where:
  !> \param[in] lun -> integer. I/O unit for error message output
  !>
  !<-----------------------------------------------------------
  subroutine kill_quadmesh(this, lun)
    implicit none
    ! vars
    class(quadmesh) :: this
    integer, intent(in) :: lun
    ! local vars
    integer :: res
    logical :: rc
    
    deallocate(&
         this%topol_inter,&
         stat=res)
    if (res.ne.0) rc=IOerr(lun, err_dealloc, 'kill_quadmesh', &
         'dealloc fail for quadmesh members topol,side_cnc,neigh,'// &
         'iside,plist')

    deallocate(&
         this%xcoord,&
         this%ycoord,&
         stat=res)
    if (res.ne.0) rc=IOerr(lun, err_dealloc, 'kill_quadmesh', &
         'dealloc fail for quadmesh members area_cell, peri_cell, bar_cell')

    call this%grid%kill(lun)

  end subroutine kill_quadmesh
  
  
  !>--------------------------------------------------------------
  !> basic quadmesh generator
  !<-----------------------------------------------------------------
  subroutine init_quadmesh(this,&
       lun_err,&
       ninter,&
       base_point, &
       measures )
    use Globals
    implicit none
    class(quadmesh),   intent(inout) :: this
    integer,           intent(in   ) :: lun_err
    integer,           intent(in   ) :: ninter
    real(kind=double), intent(in   ) :: base_point(2)
    real(kind=double), intent(in   ) :: measures(2)
    !local
    logical :: rc
    integer :: res
    integer :: i,j,k,m, nnode, nnodeincell, ncell
    integer, allocatable :: temp_topol(:,:)
    real(kind=double), allocatable :: temp_coord(:,:)

    this%grid%cell_type='square'
    this%grid%ambient_dimension=2
    this%grid%logical_dimension=2

    this%base_point = base_point
    this%measures = measures
    this%ninter   = ninter

    nnodeincell = 4
    nnode       = (ninter+1)**2
    ncell       = (ninter)**2

    allocate(&
         this%xcoord(ninter+1),&
         this%ycoord(ninter+1),&
         stat=res)
    if(res .ne. 0) rc = IOerr(lun_err, err_alloc, 'main', &
         ' member coord',res)
    
    this%xcoord(1)=base_point(1)
    this%ycoord(1)=base_point(2)
    do i=1,ninter
       this%xcoord(i+1)= base_point(1) + i* measures(1)/ninter
    end do
     do i=1,ninter
       this%ycoord(i+1)= base_point(2) + i* measures(2)/ninter
    end do

     ! topol_inter
    allocate(&
         this%topol_inter(2,ncell),&
         stat=res)
    if(res .ne. 0) rc = IOerr(lun_err, err_alloc, 'basic_quadmesh', &
         ' member topol_inter',res)
    
    m=0
    do i=1,this%ninter
       do j=1,this%ninter
          m=m+1
          this%topol_inter(1,m) = i 
          this%topol_inter(2,m) = j
       end do
    end do
    allocate(temp_coord(3,nnode), temp_topol(nnodeincell,ncell))

    m=0
    do i=1,this%ninter+1
       do j=1,this%ninter+1
          m=m+1
          temp_coord(1,m)=this%xcoord(i)
          temp_coord(2,m)=this%ycoord(j)
          temp_coord(3,m)=0
       end do
    end do
    
    
    m=0    
    do j=1,this%ninter
       do i=1,this%ninter
          m=m+1
          ! from lower left in clock wise sense
          temp_topol(1,m)=i+(j-1)*(ninter+1)
          temp_topol(2,m)=i+j*(ninter+1)
          temp_topol(3,m)=i+j*(ninter+1)+1
          temp_topol(4,m)=i+(j-1)*(ninter+1) + 1
       end do
    end do

    !
    ! build grid from data
    !
    call this%grid%init_from_data(lun_err,&
         nnode,ncell,nnodeincell,'square',&
         temp_topol, temp_coord)

    call this%grid%build_size_cell(lun_err)

    deallocate(&
         temp_coord,&
         temp_topol)

  end subroutine init_quadmesh

  
  subroutine get_interval_from_cell(this,icell, ax, bx, ay, by)
    implicit none
    class(quadmesh), intent(in) :: this
    integer,     intent(in) :: icell
    real(kind=double), intent(out) :: ax, bx, ay, by
    
    ax = this%grid%coord(1, this%grid%topol(1, icell) )
    bx = this%grid%coord(1, this%grid%topol(2, icell) )
    ay = this%grid%coord(2, this%grid%topol(2, icell) )
    by = this%grid%coord(2, this%grid%topol(3, icell) )
  end subroutine get_interval_from_cell
    

end module QuadGrid
