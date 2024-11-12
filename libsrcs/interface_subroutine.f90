subroutine refine(&
     lun_err,&
     nnodeincell,nnode,ncell, &
     coord,topol,&
     nnodeout,ncellout, &
     coordout,topolout)
  use KindDeclaration
  use Globals
  use AbstractGeometry
  use SparseMatrix
  implicit none
  integer,           intent(in   ) :: lun_err
  integer,           intent(in   ) :: nnodeincell,nnode,ncell
  integer,           intent(in   ) :: topol(nnodeincell,ncell)
  real(kind=double), intent(in   ) :: coord(3,nnode)
  integer,           intent(inout) :: nnodeout,ncellout
  integer,           intent(inout) :: topolout(nnodeincell,4*ncell)
  real(kind=double), intent(inout) :: coordout(3,3*nnode)
  ! local
  type(abs_simplex_mesh) :: grid,subgrid
  type(file) :: fgrid,fparent
  logical :: rc
  integer :: res
  type(spmat) :: connection_matrix
  character(len=256) :: cell_type
  integer, allocatable :: perm(:),iperm(:)

  if (nnodeincell==3) then
     cell_type='triangle'
  else
     cell_type='tetrahedron'
  end if

  !
  ! init grid
  !
  call grid%init_from_data(lun_err,&
       nnode,ncell,nnodeincell,cell_type, &
       topol,coord)
  write(*,*)'grid from data'
  call grid%build_size_cell(lun_err)
  write(*,*)'grid sizecell'
  call grid%build_normal_cell(lun_err)
  call grid%build_bar_cell(lun_err)  
  write(*,*)'grid done'

  write(*,*) grid%nnode
  call fgrid%init(lun_err,'gridin.dat',14,'out')
  call grid%write_mesh(lun_err,fgrid)
  call fgrid%kill(lun_err)


  !
  ! build subgrid
  !
  call subgrid%refine(lun_err,input_mesh=grid)
  write(*,*)'subgrid refined'
  call subgrid%build_nodenode_matrix(lun_err, .False.,connection_matrix)
  write(*,*)'subgrid nodenode'
  allocate(perm(subgrid%nnode),iperm(subgrid%nnode),stat=res)
  if (res.ne.0) rc = IOerr(lun_err, err_inp, 'data2grids', &
       ' work arrays perm iperm', res)
  call fparent%init(lun_err,'parent.dat',14,'out')
  call subgrid%write_parent(lun_err,fparent)
  call fparent%kill(lun_err)

  
  call connection_matrix%genrcm(6,perm,iperm)
  write(*,*)'subgrid genrcm'
  call connection_matrix%kill(lun_err)
  call subgrid%renumber(lun_err, subgrid%nnode,perm,iperm)
  write(*,*)'subgrid renumber'
  call fparent%init(lun_err,'parent_permuted.dat',14,'out')
  call subgrid%write_parent(lun_err,fparent)
  call fparent%kill(lun_err)


  call subgrid%build_size_cell(lun_err)
  call subgrid%build_normal_cell(lun_err)

  deallocate(perm,iperm,stat=res)
  if (res.ne.0) rc = IOerr(lun_err, err_dealloc, 'data2grids', &
       ' work arrays perm iperm', res)

  nnodeout=subgrid%nnode
  ncellout=subgrid%ncell

  coordout(1:3,1:nnodeout)=subgrid%coord(1:3,1:nnodeout)
  topolout(1:3,1:ncellout)=subgrid%topol(1:3,1:ncellout)
  
end subroutine refine

subroutine get_boundary_nodes_from_topol(&
     lun_err,&
     nnodeincell,nnode,ncell, &
     topol,coord,&
     nboundarynode,info,boundarynodes)
  use KindDeclaration
  use Globals
  use AbstractGeometry
  use SparseMatrix
  implicit none
  integer, intent(in   ) :: lun_err
  integer, intent(in   ) :: nnodeincell
  integer, intent(in   ) :: nnode
  integer, intent(in   ) :: ncell
  integer, intent(in   ) :: topol(nnodeincell,ncell)
  real(kind=double), intent(in   ) :: coord(3,nnode)
  integer, intent(inout) :: nboundarynode
  integer, intent(inout) :: info
  integer, intent(inout) :: boundarynodes(nnode)
  ! local
  type(abs_simplex_mesh) :: grid
  character(len=256) :: cell_type

    
   if (nnodeincell==3) then
     cell_type='triangle'
  else
     cell_type='tetrahedron'
  end if
  write(*,*) nboundarynode
  

  ! init grid
  call grid%init_from_data(lun_err,&
       nnode,ncell,nnodeincell,cell_type, &
       topol,coord)
  ! get boundary nodes
  call get_boundary_nodes(&
     lun_err,&
     grid,nboundarynode,boundarynodes)
  ! free memory
  call grid%kill(lun_err)

  write(*,*) nboundarynode,ncell

  info=0
end subroutine get_boundary_nodes_from_topol

subroutine get_boundary_nodes(&
     lun_err,&
     grid,nboundarynode,boundarynodes)
  use KindDeclaration
  use Globals
  use AbstractGeometry

  integer,                intent(in   ) :: lun_err
  type(abs_simplex_mesh), intent(inout) :: grid
  integer,                intent(inout) :: nboundarynode
  integer,                intent(inout) :: boundarynodes(grid%nnode)

  call grid%build_edge_connection(lun_err)
  call grid%build_nodebc(lun_err)
  nboundarynode=grid%nnode_bc
  do i=1,nboundarynode
     boundarynodes(i)=grid%node_bc(i)
  end do
  
end subroutine get_boundary_nodes

subroutine surrounding(&
     lun_err,&
     grid,point_coordinates, radius, &
     nsurroundingcells,surroundingcells,&
     nboundarynode,boundarynodes)
  
  use KindDeclaration
  use Globals
  use AbstractGeometry
  implicit none
  integer,                intent(in   ) :: lun_err
  type(abs_simplex_mesh), intent(inout) :: grid
  real(kind=double),      intent(in   ) :: point_coordinates(3)
  real(kind=double),      intent(in   ) :: radius
  integer,                intent(inout) :: nsurroundingcells
  integer,                intent(inout) :: surroundingcells(grid%ncell)
  integer,                intent(inout) :: nboundarynode
  integer,                intent(inout) :: boundarynodes(grid%nnode)
  ! local
  real(kind=double) :: dnrm2 
  integer :: i
  type(abs_simplex_mesh) :: grid_hole
  logical, allocatable :: selector(:)
  integer, allocatable :: map_nodes_new2old(:)

  call grid%info(6)
  write(*,*) ' here'
  call grid%build_bar_cell(lun_err)  
  
  allocate (selector(grid%ncell),map_nodes_new2old(grid%nnode))
  
  selector=.False.
  nsurroundingcells=0
  do i=1,grid%ncell
     if ( dnrm2(3,point_coordinates-grid%bar_cell(:,i),1) < radius) then
        selector(i)=.True.
        nsurroundingcells = nsurroundingcells + 1
        surroundingcells(nsurroundingcells) = i
     end if
  end do
  write(*,*) ' here',nsurroundingcells

  call grid_hole%init_selection(0, grid, grid%ncell, selector,&
       map_nodes_new2old=map_nodes_new2old)
  write(*,*) 'selected'
  call grid_hole%info(6)
  call grid_hole%build_edge_connection(lun_err)
  call grid_hole%build_nodebc(lun_err)
  write(*,*) 'nnodebc ',grid_hole%nnode_bc
  nboundarynode = grid_hole%nnode_bc
  do i = 1,nboundarynode
     write(*,*) i
     write(*,*) grid_hole%node_bc(i),  boundarynodes(i)
     boundarynodes(i)=map_nodes_new2old(grid_hole%node_bc(i))
  end do

  

  deallocate (selector)
  
end subroutine surrounding




subroutine build_refinement(&
     lun_err,&
     grid,subgrid)
  use KindDeclaration
  use Globals
  use AbstractGeometry
  use SparseMatrix
  implicit none
  integer,                intent(in   ) :: lun_err
  type(abs_simplex_mesh), intent(in   ) :: grid
  type(abs_simplex_mesh), intent(inout) :: subgrid
  ! local
  type(file) :: fgrid
  logical :: rc
  integer :: res
  type(spmat) :: connection_matrix
  character(len=256) :: cell_type
  integer, allocatable :: perm(:),iperm(:)

  !
  ! build subgrid
  !
  call subgrid%refine(lun_err,input_mesh=grid)
  write(*,*)'subgrid refined'
  call subgrid%build_nodenode_matrix(lun_err, .False.,connection_matrix)
  write(*,*)'subgrid nodenode'
  allocate(perm(subgrid%nnode),iperm(subgrid%nnode),stat=res)
  if (res.ne.0) rc = IOerr(lun_err, err_inp, 'data2grids', &
       ' work arrays perm iperm', res)
  call connection_matrix%genrcm(6,perm,iperm)
  write(*,*)'subgrid genrcm'
  call connection_matrix%kill(lun_err)
  call subgrid%renumber(lun_err, subgrid%nnode,perm,iperm)

  deallocate(perm,iperm,stat=res)
  if (res.ne.0) rc = IOerr(lun_err, err_dealloc, 'data2grids', &
       ' work arrays perm iperm', res)

end subroutine build_refinement





subroutine forcing2rhs(grid,forcing,rhs)
  use KindDeclaration
  use Globals
  use AbstractGeometry
  implicit none
  type(abs_simplex_mesh), intent(in   ) :: grid
  real(kind=double),      intent(in   ) :: forcing(grid%ncell)
  real(kind=double),      intent(inout) :: rhs(grid%nnode)
  ! local
  integer :: icell, iloc, inode
  real(kind=double) ::factorsize

  factorsize=one/grid%nnodeincell
  
  rhs = zero
  do icell = 1, grid%ncell
     do iloc = 1,grid%nnodeincell
        
        inode = grid%topol(iloc,icell)
        rhs(inode) =  rhs(inode) + &
             factorsize * forcing(icell) * grid%size_cell(icell)
     end do
  end do
  
end subroutine forcing2rhs


subroutine build_subgrid_rhs(subgrid,rhs,forcing, dirac)
  use KindDeclaration
  use Globals
  use AbstractGeometry
  implicit none
  type(abs_simplex_mesh),      intent(in   ) :: subgrid
  real(kind=double),           intent(inout) :: rhs(subgrid%nnode)
  real(kind=double),  intent(in   ) :: forcing(subgrid%ncell_parent)
  real(kind=double),  intent(in   ) :: dirac(subgrid%nnode_parent)
  
  ! local
  logical :: rc
  integer :: res,n1,n2,inode
  real(kind=double), allocatable :: forcing_subgrid(:)
  real(kind=double), allocatable :: dirac_subgrid(:)

  allocate(forcing_subgrid(subgrid%ncell),dirac_subgrid(subgrid%nnode),stat=res)
  if (res.ne.0) rc = IOerr(0, err_alloc, 'build_subgrid_rhs', &
       ' work arrays forcing_subgrid,dirac_subgrid', res)

  rhs=zero
  !
  ! interpolate data to subgrid
  !
  write(*,*) 'forcing'
  forcing_subgrid=zero

  call subgrid%proj_subgrid( forcing, forcing_subgrid)
  write(*,*) 'proj'
  call forcing2rhs(subgrid,forcing_subgrid,rhs)
  write(*,*) 'rhs'
  
  dirac_subgrid=zero
  do inode=1,subgrid%nnode
     n1=subgrid%node_parent(1,inode)
     n2=subgrid%node_parent(2,inode)
     if (n1 .eq. n2 ) then
        dirac_subgrid(inode)=dirac(n1)
     end if
  end do
  write(*,*) 'dirac'


  !
  ! buils rhs 
  !
  rhs=rhs+dirac_subgrid
  
  deallocate(forcing_subgrid,dirac_subgrid,stat=res)
  if (res.ne.0) rc = IOerr(0, err_dealloc, 'build_subgrid_rhs', &
       ' work arrays forcing_subgrid,dirac_subgrid', res)

  
end subroutine build_subgrid_rhs


! subroutine uniform_refinement(lun_err,nref,grid,subgrid)
!   use Globals
!   use AbstractGeometry
!   use SparseMatrix

!   implicit none
!   integer,           intent(in   ) :: lun_err
!   integer,           intent(in   ) :: nref
!   integer,           intent(in   ) :: nnodeincell,nnode,ncell
!   integer,           intent(in   ) :: topol(nnodeincell,ncell)
!   real(kind=double), intent(in   ) :: coord(3,nnode)
!   integer,           intent(inout) :: nnodeout,ncellout
!   integer,           intent(inout) :: topolout(nnodeincell,4*ncell)
!   real(kind=double), intent(inout) :: coordout(3,3*nnode)
!   type(file) :: fgrid
!   type(file) :: fsubgrid
!   type(file) :: fparent
!   integer  :: nref
!   logical  :: rc
!   integer  :: res
!   character(len=256) :: input,path_grid,path_subgrid,path_parent

!   integer :: stderr,stdout,debug
!   integer :: i, itria,err,narg,nargs
!   integer :: ngrids,flag_reorder

!   type(spmat) :: adj_matrix
!   integer, allocatable :: perm(:),iperm(:)

!   stderr=0
!   stdout=6
!   debug=6

!   if (nnodeincell==3) then
!      cell_type='triangle'
!   else
!      cell_type='tetrahedron'
!   end if
  
!   allocate(grids(0:nref),stat=res)
!   if(res .ne. 0) rc = IOerr(stderr, err_alloc, 'main', &
!        ' subgrids',res)

!   call grid(0)%init_from_data(lun_err,&
!        nnode,ncell,nnodeincell,cell_type, &
!        topol,coord)

!   write(stdout,'(a,I2,a,I8,a,I8)') ' Refinement # = ', 0,&
!        ' Nnode= ', grids(0)%nnode,&
!        ' Ncell= ', grids(0)%ncell
  
!   !
!   ! refineing procedure
!   !    
!   do i=1,nref 
!      call grids(i-1)%build_edge_connection(stderr)
!      call grids(i)%refine(stderr,input_mesh=grids(i-1))
!      write(stdout,'(a,I2,a,I8,a,I8)') ' Refinement # = ', i,&
!           ' Nnode= ', grids(i)%nnode,&
!           ' Ncell= ', grids(i)%ncell
!   end do

!   select case (flag_reorder)
!   case (1)
!      !
!      ! build adjency matrix
!      !
!      call grids(nref)%build_nodenode_matrix(stderr,.False.,adj_matrix)
!      adj_matrix%coeff=one

!      !
!      ! build rcm permutation of graph node-node
!      !
!      allocate(perm(grids(nref)%nnode),iperm(grids(nref)%nnode))
!      call adj_matrix%genrcm(6,perm,iperm)

!      !
!      ! reorder
!      !
!      call grids(nref)%renumber(stderr, grids(nref)%nnode,perm,iperm)    

!      !
!      ! free memory
!      !
!      call adj_matrix%kill(stderr)
!      deallocate(perm,iperm,stat=res)

!   case (2)
!      if ( grids(nref)%nnodeincell .eq. 3 ) then
!         call grids(nref)%build_cellcell_matrix(stderr,.False.,adj_matrix)
!         adj_matrix%coeff=one

!         !
!         ! build rcm permutation of graph node-node
!         !
!         allocate(perm(grids(nref)%ncell),iperm(grids(nref)%ncell))
!         call adj_matrix%genrcm(6,perm,iperm)

!         !
!         ! reorder
!         !
!         call grids(nref)%renumber(stderr, grids(nref)%ncell,perm,iperm)    

!         !
!         ! free memory
!         !
!         call adj_matrix%kill(stderr)
!      else
!         rc = IOerr(stderr, err_inp, 'main', &
!              ' Cell Renumbering defined only for triangles =')
!      end if
!   case (3)
!      !
!      ! build adjency matrix
!      !
!      call grids(nref)%build_nodenode_matrix(stderr,.False.,adj_matrix)
!      adj_matrix%coeff=one

!      !
!      ! build rcm permutation of graph node-node
!      !
!      allocate(perm(grids(nref)%nnode),iperm(grids(nref)%nnode))
!      call adj_matrix%genrcm(6,perm,iperm)

!      !
!      ! reorder
!      !
!      call grids(nref)%renumber(stderr, grids(nref)%nnode,perm,iperm)    

!      !
!      ! free memory
!      !
!      call adj_matrix%kill(stderr)
!      deallocate(perm,iperm,stat=res)

!      !
!      ! cell reorder
!      !
!      if ( grids(nref)%nnodeincell .eq. 3 ) then
!         call grids(nref)%build_cellcell_matrix(stderr,.False.,adj_matrix)
!         adj_matrix%coeff=one

!         !
!         ! build rcm permutation of graph node-node
!         !
!         allocate(perm(grids(nref)%ncell),iperm(grids(nref)%ncell))
!         call adj_matrix%genrcm(6,perm,iperm)

!         !
!         ! reorder
!         !
!         call grids(nref)%renumber(stderr, grids(nref)%ncell,perm,iperm)    

!         !
!         ! free memory
!         !
!         call adj_matrix%kill(stderr)
!      else
!         rc = IOerr(stderr, err_inp, 'main', &
!              ' Cell Renumbering defined only for triangles =')
!      end if
!   end select

!   nnodeout=subgrid%nnode
!   ncellout=subgrid%ncell

!   coordout(1:3,1:nnodeout)=subgrid%coord(1:3,1:nnodeout)
!   topolout(1:3,1:ncellout)=subgrid%topol(1:3,1:ncellout)

!   !
!   ! free memory
!   !
!   do i=1,nref 
!      call grids(i)%kill(stderr)
!   end do  


! end subroutine uniform_refinement
