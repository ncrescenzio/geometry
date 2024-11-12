!>-------------------------------------------------------------
!> defines objects containing geometry and mesh information
!>
!> objects:
!> type mesh --> TODO Description
!>
!> methods:
!> 
!<------------------------------------------------------------- 
module Geometry3d
  use Globals
  use AbstractGeometry
  implicit none
  private
  !> structure containing variable mesh and geometrical information 
  type, extends(abs_simplex_mesh),public :: mesh3d
     !> Number of mesh nodes.
     !> **Input in the constructor**.
     !integer :: nnode = 0  
     !> Number of mesh tetrahedrons.
     !> **Input in the constructor**.
     !integer :: ncell = 0 
     !> Number of nodes in each cell.
     !> **Input in the constructor**.
     !integer  :: nnodeincell = 4 
     !> Number of mesh faces (including boundary). 
     !> Calculated in the constructor.
     !integer :: nface = 0
     !> Number of boundary faces.
     !> **Calculated in the constructor.**
     !integer :: nface_bc = 0
     !> Number of mesh edges (including boundary). 
     !> Calculated in the constructor.
     !integer :: nedge = 0 
     !> Number of boundary edges.
     !> **Calculated in the constructor.**
     !integer :: nedge_bc = 0
     !> Number of boundary node.
     !> **Calculated in the constructor.**
     !integer :: nnode_bc = 0 
     !> level of mesh topological quantities optionally calculated
     !> (true/false)
     !> (default =0, =1 first level cnc, =2 second level cnc)
     !>
     !> Default level cnc vars:
     !> topol, coord, size_cell, bar_cell
     !> First level cnc vars:
     !> face_plist, face_cnc, edge_cnc, face_iside, iside_edge, neigh_face,
     !> surface_tetra, length_edge, bar_edge, normal
     integer :: cnc_built = 0
     !>-------------------------------------------------------------------------
     !> Topological info
     !> Dimension (5,ncell).
     !> **Input in the constructor.**
     !integer, allocatable :: topol(:,:)
     !> Dimension (4,ncell). 
     !> Faces (in the global enumeration) in each cell.
     !> **Calculated in the constructor.**
     integer, allocatable :: side_cnc(:,:)
     !> Dimension (6,ncell). 
     !> Edges (in the global enumeration) in each cell.
     !> **Calculated in the constructor.**
     !integer, allocatable :: edge_cnc(:,:)
     !> Dimension (4,ncell).
     !> Neighboring tetrahedrons for each tetrahedron.
     !> **Calculated in the constructor.**
     !integer, allocatable :: neigh_face(:,:)
     !> Dimension (3,nface).
     !> Nodes for each face. 
     !> **Calculated in the constructor.**
     integer, allocatable :: iside_face(:,:)
     !> Dimension (2,nedge). 
     !> Nodes for each edge. 
     !> **Calculated in the constructor.**
     integer, allocatable :: iside_edge(:,:)
     !> Dimension (nnode_bc). 
     !> Nodes at boundary 
     !> **Calculated in the constructor.**
     !integer, allocatable :: node_bc(:)
     !> Dimension (2,nface). 
     !> East and West cell for each face
     !> East and West are calculated with respect to the nodal ordering in CELL.
     !> **Calculated in the constructor.**
     integer, allocatable :: plist(:,:)
     !> Dimension (3,nnode). 
     !> Nodal coordinate. 
     !> **Input in the constructor.**
     !real(kind=double), allocatable :: coord(:,:)
     !> Dimension (ncell).
     !> Volume of the cell
     !> **Computed in procedure Prop.**
     !real(kind=double), allocatable :: size_cell(:)
     !> Dimension (ncell).
     !> Surface of the cell
     !> **Computed in procedure Prop.**
     real(kind=double), allocatable :: surface_tetra(:)
     !> Dimension (3,ncell).
     !> Coordinates of center of mass of the cell
     !> **Computed in procedure Prop.**
     !real(kind=double), allocatable :: bar_cell(:,:)
     !> Dimension (nface).
     !> Surface of the face
     !> **Computed in procedure Prop.**
     real(kind=double), allocatable :: surface_face(:)
     !> Dimension (3,nface).
     !> Normal vector to the surface (consistent with face_iside and face_plist)
     !> **Computed in procedure Prop.**
     real(kind=double), allocatable :: normal(:,:)
     !> Dimension (3,nedge). 
     !> Coord of gravity center of edges
     !> **Computed in procedure Prop.**
     real(kind=double), allocatable :: bar_edge(:,:)
     !> Dimension (nnode)
     !> permutation vector for nodes
     !> **Calculated in the constructor.**
     integer, allocatable :: perm_nodes(:)
     !>-------------------------------------------------------------------------
     !> Arrays allocated only for grid levels (using uniform refinements) 
     !> (grid level 0 is the coarsest)
     !>-------------------------------------------------------------------------
     !> Grid level 
     !integer  :: grid_level=0
     !> Number of nodes of parent grid
     !integer :: nnode_parent=0
     !> Number of tetrahedrons parent grid = ncell/8
     !integer :: ncell_parent=0
!!$     !> Number of edges of parent grid
!!$     integer :: nedge_parent=0
     !> Parent nodes of each node (can be the same)
     !> Dimension (2,nnode)
     !> A node can be generated by a single parent
     !> (if it corresponds to a vertex node of the parent mesh)
     !> or by two parents (if it is an edge midpoint node of the
     !> parent mesh).
     !> E.G.: if (node_parent(1,inode) .eq. node_parent(2,inode) )
     !> inode has only one parent with index node_parent(1,inode)
     !> (is its direct son),
     !> otherwise inode is the mid point of the edge defined by
     !> node_parent(1,inode) and node_parent(2,inode)
     !> (this type of array and of reasoning is needed for nodes and
     !> edges) 
     !integer, allocatable :: node_parent(:,:)
     !> Dimension (ncell)
     !> Parent tetrahedrons
     !integer, allocatable :: cell_parent(:)
   contains
     !> static constructor
     !> (procedure public for type mesh3d)
     procedure, public, pass :: init => init_mesh
     !> static destructor
     !> (procedure public for type mesh3d)
     !procedure, public, pass :: kill => kill_mesh
     !> Info procedure.
     !> (public procedure for type mesh3d)
     procedure, public, pass :: info => info_mesh
     !> Procedure to write mesh3d in read-ready format
     !> (public procedure for type mesh3d)
     procedure, public, pass :: write => write_mesh
     !> Reads mesh3d from files
     !> (public procedure for type mesh3d)
     procedure, public, pass :: read_mesh_3d
     !> builds mesh3d refining a previous mesh 
     !> (private procedure for type mesh)
     !procedure, public, pass :: refine => build_refined_topol_coord
     !> Build mesh faces connectivity data starting from topol
     !> (private procedure for type mesh, used in init)
     !procedure, public, pass :: connection_faces
     !> Build mesh edges connectivity data starting from topol
     !> (private procedure for type mesh, used in init)
     procedure, public, pass :: connection_edges
     !> Calculate cell areas, edge lengths, normals, etc
     !> (private procedure for type mesh, used in init)
     procedure, private, pass :: meshProp
     !> calculate the average of real array defined over subgrid
     !> (public procedure for type mesh)
     procedure, public, pass :: avg
     !> calculate the average of real array defined over subgrid
     !> (public procedure for type mesh)
     procedure, public, pass :: avg_vec
     !> calculate the average of real array defined over subgrid
     !> (public procedure for type mesh)
     procedure, public, pass :: proj_subgrid
     !> Project data defined on nodes to subgrid nodes
     !> (public procedure for type mesh)
     procedure, public, pass :: projnode_subgrid
     !> calculate ancestrors in series of refined meshes
     !> (public procedure for type mesh)
     procedure, public, pass :: heraldry
     !> calculate Lebesgue p-norm of piecewise-constant array
     !> (public procedure for type mesh)
     procedure, public, pass :: normp_cell
     !> calculate parameter of the spatial discretization
     !> (public procedure for type mesh)
     procedure, public, pass :: meshpar
     !> procedure to find the index of the node
     !> closer to a given point
     procedure, public, pass :: closer_node
     !> procedure to find the index of the node
     !> closer to a given point
     procedure, public, pass :: info_point
     !> Renumber node numbering for 
     !> (public procedure for type mesh)
     procedure, public, pass :: reorder
     !> Renumber node numbering for 
     !> (private procedure for type mesh)
     procedure, private, pass :: build_rcm_permutation
     !> calculate the projection along the edge outward unit normal
     !> (public procedure for type mesh)
     procedure, public, pass :: proj => normal_projection
     !> Procedure to read subgrid cell_parent and node_parent
     !> (public procedure for type mesh)
     procedure, public, pass :: read_parent => read_grid_subgrid_relation
     !> Procedure to write subgrid cell_parent and node_parent
     !> (public procedure for type mesh)
     procedure, public, pass :: write_parent => write_grid_subgrid_relation
     !> Procedure to check grid subgrid dependence
     !> (public procedure for type mesh)
     procedure, public, pass :: check_relations => check_subgrid_grid_relation
     !> Procedure to compute node2node_connection
     !> (public procedure for type mesh)
     procedure, public, pass :: nodenode_connection
     !> Procedure for building
     !> $\int _{Cell_r} |\Laplacain \Pot | 
     !> (procedure public for type tdpotsys)
     procedure, public , pass :: p0_laplacian
  end type mesh3d
contains
  !>-------------------------------------------------------------
  !> Static constructor.
  !> (procedure public for type mesh)
  !> Instantiate (allocate if necessary)
  !> and initialize (by also reading from input file if requested)
  !> variable of type mesh
  !>
  !> usage:
  !>     call 'var'%init(lun_err,&
  !>                flag_read_build, &
  !>                [file2read[=...], &
  !>                input_mesh[=...], &
  !>                flag_cnc[=...], &
  !>                flag_reorder[=...]
  !>                ])
  !>
  !> where:
  !> \param[in] lun_err                 -> integer. IO unit for error msg
  !> \param[in] flag_read_build         -> integer. flag=0 reads from file
  !>                                                flag=1 refine from input_mesh
  !>                                                flag=2 read refined mesh from file
  !> \param[in] (optional) grid2read    -> type(file). I/O file for mesh data
  !> \param[in] (optional) parent2read  -> type(file). I/O file for parent data
  !> \param[in] (optional) input_mesh   -> type(mesh3d). mesh to refine
  !>                                       (used only if flag_read_build=1)
  !> \param[in] (optional) flag_cnc     -> integer.
  !>                                       flag=0 build only topol, coord,
  !>                                              and related real properties
  !>                                       flag=1 in addition, build
  !>                                              level 1 connectivity vars, 
  !>                                              and related real properties
  !>                                       flag=2 in addition, build
  !>                                              level 2 connectivity vars, 
  !>                                              and related real properties
  !> \param[in] (optional) flag_reorder -> integer.
  !>                                       flag=0 build mesh only
  !>                                       flag=1 build mesh and do reorder
  !>                                              minimizing the bandwith
  !>                                              (Cathill-McKee) for P1 Galerkin FEM
  !<-------------------------------------------------------------
  subroutine init_mesh(this,lun_err,&
       flag_read_build, & ! next are optional arguments
       grid2read, parent2read, input_mesh,&
       flag_cnc, flag_reorder)
    use Globals
    implicit none
    !vars
    class(mesh3d),           intent(inout) :: this
    integer,               intent(in   ) :: lun_err
    integer,               intent(in   ) :: flag_read_build
    type(file),  optional, intent(in   ) :: grid2read, parent2read
    type(mesh3d),  optional, intent(in   ) :: input_mesh
    integer,     optional, intent(in   ) :: flag_cnc, flag_reorder

    !local 
    logical :: rc
    integer :: res, tmp_flag
    integer, allocatable :: perm(:), inv_perm(:)
    integer :: i,j

    this%nnodeincell = 4
    this%ambient_dimension = 3
    this%logical_dimension = 3
    this%cell_type = 'tetrahedron'
    this%cell_id = 10
    this%mesh_type = '3d'
    
    ! Initialized topol and coord
    select case (flag_read_build)
    case (0)
       ! Read coord and topol from file
       call this%read_mesh_3d(lun_err,grid2read)
    case (1)
       ! build uniformly refined mesh
       call this%refine(lun_err,input_mesh)
    case (2)
       ! read refined mesh+parent
       call this%read_mesh_3d(lun_err,grid2read)
       call this%read_parent(lun_err,parent2read)
    case (3) 
       ! do nothing nnode, ncell, topol and coord are already defined
       ! move to the connections part
    case default
       write(lun_err,*) ' Wrong Flag in procedure init mesh'
       write(lun_err,*) ' Flag_read_build not equal to 0 or 1 or 2'
       stop
    end select
    ! set the number of zones
    this%nzone = maxval(this%topol(this%nnodeincell+1,:),this%ncell)

    if ( present(flag_reorder) .and. (flag_reorder .eq. 1)  ) then
       ! allocation of local working arrays
       allocate(perm(this%nnode),inv_perm(this%nnode),stat=res)
       if(res .ne. 0) rc = IOerr(lun_err,err_alloc ,'init_mesh',&
            'work array perm, inv_perm',res)
       ! Compute node permutation 
       call this%build_rcm_permutation(lun_err,perm,inv_perm)
       ! Permute nodes
       call this%reorder(lun_err,perm,inv_perm)
       ! Create a copy 
       allocate(this%perm_nodes(this%nnode),stat=res)
       if(res .ne. 0) rc = IOerr(lun_err,err_alloc,'init_mesh',&
            ' perm_nodes',res)
       this%perm_nodes = perm
       ! Free memory
       deallocate(perm,inv_perm,stat=res) 
       if(res .ne. 0) rc = IOerr(lun_err,err_dealloc,'init_mesh',&
            'work array perm',res) 
    end if
    
    ! Build connection_faces and connection_edges 
    if ( present(flag_cnc) ) this%cnc_built=flag_cnc
    select case (this%cnc_built)
    case default
       call this%meshProp(lun_err)
    case (1)
       call this%build_face_connection(lun_err)
       call this%connection_edges(lun_err)
       call this%meshProp(lun_err)
    end select

    do i=1,this%ncell
       call isort(4,this%topol(1:4,i))
    end do

  end subroutine init_mesh
  !>-------------------------------------------------------------
  !> Part of the static constructor.
  !> (procedure private for type mesh)
  !> Instantiate (allocate if necessary)
  !> and initialize (by also reading from input file)
  !> variable of type mesh
  !>
  !> usage:
  !>     call this%read_mesh(lun_err, file2read)
  !>
  !> where:
  !> \param[in] lun_err   -> integer. IO unit for error msg
  !> \param[in] file2read -> type(file). I/O file information
  !<-------------------------------------------------------------
  subroutine read_mesh_3d(this, lun_err, file2read)
    use Globals
    implicit none
    !vars
    class(mesh3d),    intent(inout) :: this
    integer,        intent(in ) :: lun_err
    type(file),     intent(in ) :: file2read
    ! local vars
    integer :: u_number
    integer :: j, k, res, nnodeincell    
    logical :: rc
    character(len=256) :: rdwr,str,fname

    u_number  = file2read%lun
    fname     = file2read%fn
    nnodeincell = 4

    read(u_number,*,iostat=res) this%nnode
    if(res .ne. 0) rc = IOerr(lun_err, err_inp , 'read_mesh', &
         etb(fname) // ' vars nnode',res)

    read(u_number,*,iostat=res) this%ncell
    if(res .ne. 0) rc = IOerr(lun_err, err_inp , 'read_mesh', &
         etb(fname) // ' vars ncell',res)

    allocate(this%coord(3,this%nnode),stat=res)
    if(res .ne. 0) rc = IOerr(lun_err, err_alloc, 'read_mesh', &
         '  type mesh member coord (array)',res)
    allocate(this%topol(5,this%ncell),stat=res)
    if(res .ne. 0) rc = IOerr(lun_err, err_alloc, 'read_mesh', &
         '  type mesh member tetrahedrons (array)',res)

    do j=1,this%nnode
       read(u_number,*,iostat=res) (this%coord(k,j),k=1,3)
       if(res .ne. 0) THEN
          write(rdwr,'(i5)') j
          str=etb(rdwr)//'/'
          rc = IOerr(lun_err, err_inp , 'read_mesh', &
               trim(etb(fname)) //&
               ' type mesh member array coord at line '//&
               trim(str),res)
       end if
    end do
    

    do j=1,this%ncell
       read(u_number,*,iostat=res) (this%topol(k,j), k=1,this%nnodeincell+1)
       if(res .ne. 0) then
          write(rdwr,'(i5)') j
          str=etb(rdwr)//'/'
          write(rdwr,'(i5)') k
          str=trim(str)//etb(rdwr)
          rc = IOerr(lun_err, err_inp , 'read_mesh', &
               trim(etb(fname)) // &
               ' type mesh member array tetra  at line/col '//trim(str),res)
       end if
    end do
    this%nzone = maxval(this%topol(nnodeincell+1,:),this%ncell)

  end subroutine read_mesh_3d
  !>-------------------------------------------------------------
  !> Procedure to read parent grid information from file
  !> (procedure public for type mesh)
  !> Instantiate (allocate if necessary)
  !> and initialize (by also reading from input file)
  !> variable of type mesh
  !>
  !> usage:
  !>     call 'var'%read_parent(lun_err, file2read)
  !>
  !> where:
  !> \param[in] lun_err    -> integer. IO unit for error msg
  !> \param[in] file2read  -> type(file). I/O file information
  !>
  !> Remark: the first number read in is this%grid_level. This refers
  !> to the grid level of the parent mesh contained in file2read.
  !<-------------------------------------------------------------
  subroutine read_grid_subgrid_relation(this, lun_err, file2read)
    use Globals
    implicit none
    !vars
    class(mesh3d),    intent(inout) :: this
    integer,        intent(in   ) :: lun_err
    type(file),     intent(in   ) :: file2read

    !local 
    logical :: rc
    integer :: res
    integer :: k, inode, icell, parent_level
    integer :: u_number
    character(len=256) :: rdwr,str,fname
    integer :: nnode_parent_read,ncell_parent_read
    integer :: subnnode_read,subncell_read

    u_number = file2read%lun
    fname = file2read%fn

    read(u_number,*,iostat=res) parent_level
    if(res .ne. 0) rc = IOerr(lun_err, err_inp , 'read_grid_subgrid_relation', &
         etb(fname) // ' vars nnode',res)
    this%grid_level = parent_level + 1
    
    ! read  array dimension and check that the mesh is consistent
    read(u_number,*,iostat=res) nnode_parent_read,  subnnode_read 
    if(res .ne. 0) rc = IOerr(lun_err, err_inp , 'read_grid_subgrid_relation', &
         etb(fname) // ' vars nnode_parent_read,  subnnode_read ',res)
    if ( subnnode_read .ne. this%nnode) then
       write(lun_err,*) ' Mismatch node number'
       write(lun_err,*) ' subgrid nnode = ', this%nnode
       write(lun_err,*) ' subnnode read  = ', subnnode_read 
       stop
    end if

    read(u_number,*,iostat=res) ncell_parent_read,  subncell_read 
    if(res .ne. 0) rc = IOerr(lun_err, err_inp , 'read_grid_subgrid_relation', &
         etb(fname) // ' vars ncell_parent_read, subncell_read',res)
    if ( subncell_read .ne. this%ncell) then
       write(lun_err,*) ' Mismatch node number'
       write(lun_err,*) ' subgrid ncell = ', this%ncell
       write(lun_err,*) ' subncell read = ', subncell_read 
       stop
    end if

    ! set dimension arrays
    this%nnode_parent = nnode_parent_read
    this%ncell_parent = ncell_parent_read
    ! set allocation
    allocate(this%node_parent(2,this%nnode),stat=res)
    if(res .ne. 0) rc = IOerr(lun_err, err_alloc, 'read_grid_subgrid_relation', &
         '  type mesh member node_parent',res)
    allocate(this%cell_parent(this%ncell),stat=res)
    if(res .ne. 0) rc = IOerr(lun_err, err_alloc, 'read_grid_subgrid_relation', &
         '  type mesh member cell_parent',res)

    do inode=1,this%nnode
       read(u_number,*,iostat=res) (this%node_parent(k,inode),k=1,2)
       if(res .ne. 0) then
          write(rdwr,'(i5)') inode
          str=etb(rdwr)//'/'
          rc = IOerr(lun_err, err_inp , 'read_grid_subgrid_relation', &
               trim(etb(fname)) //&
               ' type mesh member array coord at line '//&
               trim(str),res)
       end if
    end do

    do icell=1,this%ncell
       read(u_number,*,iostat=res) this%cell_parent(icell)
       if(res .ne. 0) then
          write(rdwr,'(i5)') icell
          str=etb(rdwr)
          rc = IOerr(lun_err, err_inp , 'read_grid_subgrid_relation', &
               trim(etb(fname)) // &
               ' type mesh member array cell_parent at line/col '//trim(str),res)
       end if
    end do
  end subroutine read_grid_subgrid_relation
  !>-------------------------------------------------------------
  !> Procedure to write grid parent information to file
  !> (procedure public for type mesh)
  !> Instantiate (allocate if necessary)
  !> and initialize (by also reading from input file)
  !> variable of type mesh
  !>
  !> usage:
  !>     call 'var'%write_parent(lun_err, file2write)
  !>
  !> where:
  !> \param[in] lun_err    -> integer. IO unit for error msg
  !> \param[in] file2write -> type(file). I/O file information
  !>
  !> Remark: the first number written out is this%grid_level. This refers
  !> to the grid level of the parent mesh contained in file2write.
  !<-------------------------------------------------------------
  subroutine write_grid_subgrid_relation(this, lun_err, file2write)
    use Globals
    implicit none
    !vars
    class(mesh3d),    intent(in) :: this
    integer,        intent(in) :: lun_err
    type(file),     intent(in) :: file2write
    !local 

    integer :: k, inode, icell
    integer :: u_number
    character(len=256) :: rdwr,str,fname
    logical :: rc
    integer :: res


    u_number = file2write%lun
    fname = file2write%fn
    write(u_number,*,iostat=res) this%grid_level
    if(res .ne. 0) rc = IOerr(lun_err, err_out , 'write_grid_subgrid_relation', &
         etb(fname) // ' vars grid_level',res)

    if (this%grid_level.eq.0) then
       rc=IOerr(lun_err, wrn_out , &
            'write_grid_subgrid_relation', &
            etb(fname) // ' grid_level 0: no parent info',res)
       return
    end if

    ! Write array dimensions
    write(u_number,*,iostat=res) this%nnode_parent, this%nnode
    if(res .ne. 0) rc = IOerr(lun_err, err_out , 'write_grid_subgrid_relation', &
         etb(fname) // ' vars nnode',res)

    write(u_number,*,iostat=res) this%ncell_parent, this%ncell
    if(res .ne. 0) rc = IOerr(lun_err, err_out, 'write_grid_subgrid_relation', &
         etb(fname) // ' vars ncell',res)

    ! Write node parent
    do inode=1,this%nnode
       write(u_number,*,iostat=res) (this%node_parent(k,inode),k=1,2)
       if(res .ne. 0) then
          write(rdwr,'(i5)') inode
          str=trim(adjustl(rdwr))//'/'
          rc = IOerr(lun_err, err_out , 'write_grid_subgrid_relation', &
               trim(etb(fname)) //&
               ' type mesh member array coord at line '//&
               trim(str),res)
       end if
    end do

    ! Write tetra parent
    do icell=1,this%ncell
       write(u_number,*,iostat=res) this%cell_parent(icell)
       if(res .ne. 0) then
          write(rdwr,'(i6)') icell
          str=trim(adjustl(rdwr))
          rc = IOerr(lun_err, err_out , 'write_grid_subgrid_relation', &
               trim(etb(fname)) // &
               ' type mesh member array tetra at tetra nmb: '//trim(str),res)
       end if
    end do

  end subroutine write_grid_subgrid_relation
  !>-------------------------------------------------------------
  !> Procedure to check grid/subgrid relations
  !> (procedure public for type mesh)
  !> Return ierr not equal to zero if something is not right
  !> 
  !> usage:
  !>     ierr=check_relations(grid,subgrid,flag)
  !>
  !> where:
  !> \param[in] grid     -> type(mesh3d). 
  !> \param[in] subgrid  -> type(mesh3d). 
  !> \param[in] flag     -> integer. flag=1 check grid_level
  !>                                 flag=2 check coordinates
  !> \param[out] ierr    -> integer. ierr=1 wrong level of subgrid
  !>                                 ierr=2 wrong cordinate
  !<-------------------------------------------------------------
  function check_subgrid_grid_relation(subgrid,grid,flag) result (ierr)
    use Globals
    implicit none
    ! External variables
    class(mesh3d),  intent(in) :: subgrid
    type(mesh3d),   intent(in) :: grid
    integer,      intent(in) :: flag
    integer :: ierr
    !local 
    integer :: inode,n1,n2
    real(kind=double) :: temp(3),dnrm2
      
    ierr= 0
    if ( flag .ge. 1 )  then
       if ( (grid%grid_level .ne. subgrid%grid_level -1) .or. &
            ( subgrid%ncell/8 .ne. grid%ncell) ) then
          ierr=1
       end if
       if ( ( ierr .eq. 0) .and. (flag .eq.  2) ) then
          do inode = 1, subgrid%nnode
             n1 = subgrid%node_parent(1,inode)
             n2 = subgrid%node_parent(2,inode)
             
             temp(:) = subgrid%coord(:,inode) - &
                  onehalf*(grid%coord(:,n1) +grid%coord(:,n2))
             if ( dnrm2(3,temp,1) .gt. 1.0d-13 ) then
                ierr=2
                exit
             end if
          end do
       end if
    end if
  end function check_subgrid_grid_relation
  !>-------------------------------------------------------------
  !> Part of the build_mesh procedure
  !> (procedure private for type mesh)
  !> Build the topol and coord array of a conformally refined
  !> grid from a given coarser grid.
  !> Builds also the parent relationships.
  !>
  !> usage:
  !>     call 'var'%refine(lun_err, input_mesh)
  !>
  !> where:
  !> \param[in ] lun_err      -> integer. IO unit for err. msg. 
  !> \param[in ] input_mesh -> type(mesh3d). mesh to be refined
  !>
  !<-------------------------------------------------------------
  subroutine build_refined_topol_coord(this, lun_err, input_mesh)
    use Globals
    implicit none
    ! External variables
    class(mesh3d), intent(inout) :: this
    integer,     intent(in)  :: lun_err
    type(mesh3d),  intent(in)  :: input_mesh
    ! Local vars
    type(mesh3d) :: local_mesh 
    integer :: vert(4),edg(3),midnod(2),nod(3,2), nodes_local(4,4)
    integer :: iedge,inode,icell,ivert,izone,itetloc
    integer :: res
    integer :: iedgeloc,jedgeloc,inodloc,i,ntetref(8),iloc
    logical :: rc
    integer :: n1,n2,n3,n4,n12,n13,n14,n23,n24,n34
    real(kind=double) :: norm1,norm2,norm3
    real(kind=double) :: difference1_x,difference1_y,difference1_z
    real(kind=double) :: difference2_x,difference2_y,difference2_z
    real(kind=double) :: difference3_x,difference3_y,difference3_z
    
    local_mesh=input_mesh
    if (local_mesh%cnc_built .lt. 1) then
       call local_mesh%build_face_connection(lun_err)
       call local_mesh%connection_edges(lun_err)
       local_mesh%cnc_built = 1
       call local_mesh%meshProp(lun_err)
    end if

    ! assign subgrid dimension
    this%nnode        = local_mesh%nnode + local_mesh%nedge
    this%ncell        = 8 * local_mesh%ncell
    this%grid_level   = local_mesh%grid_level + 1
    this%nnode_parent = local_mesh%nnode
    this%ncell_parent = local_mesh%ncell

    ! allocate additional array for subgrid
    allocate (&
         this%topol(5,this%ncell),&
         this%coord(3,this%nnode),&
         this%cell_parent(this%ncell),&
         this%node_parent(2,this%nnode),&
         stat=res)
    if(res .ne. 0) rc = IOerr(lun_err, err_alloc, 'build_refined_topol_coord', &
         ' type mesh member topol, coord, cell_parent'//&
         'node_parent',res)
    
    ! assign coordinates
    this%coord(:,1:local_mesh%nnode)=local_mesh%coord
    do iedge=1,local_mesh%nedge
       this%coord(:,local_mesh%nnode+iedge)=local_mesh%bar_edge(:,iedge)
    end do

    ! subnode parent
    do inode = 1, local_mesh%nnode
       this%node_parent(1,inode) = inode  
       this%node_parent(2,inode) = inode
    end do
    do iedge = 1, local_mesh%nedge 
       inode = local_mesh%nnode + iedge
       this%node_parent(1,inode) = local_mesh%iside_edge(1,iedge)
       this%node_parent(2,inode) = local_mesh%iside_edge(2,iedge)
    end do

    ! initialize nodes_local and ntetref
    nodes_local = 0
    do icell=1,8
       ntetref(icell) = icell
    end do

    ! build new topol
    do icell=1,local_mesh%ncell
       izone = local_mesh%topol(5,icell)
       ! sort the element nodes of coarser mesh  in increasing order
       call isort(4,local_mesh%topol(1:4,icell))
       ! nodes new topol
       n1  = local_mesh%topol(1,icell)
       n2  = local_mesh%topol(2,icell)
       n3  = local_mesh%topol(3,icell)
       n4  = local_mesh%topol(4,icell)
       n12 = local_mesh%nnode + local_mesh%edge_cnc(1,icell)
       n13 = local_mesh%nnode + local_mesh%edge_cnc(2,icell)
       n14 = local_mesh%nnode + local_mesh%edge_cnc(3,icell)
       n23 = local_mesh%nnode + local_mesh%edge_cnc(4,icell)
       n24 = local_mesh%nnode + local_mesh%edge_cnc(5,icell)
       n34 = local_mesh%nnode + local_mesh%edge_cnc(6,icell)
!!$       write(*,*) n1,n2,n3,n4,n12,n13,n14,n23,n24,n34
       ! matrix local nodes
       nodes_local(1,1) = n1
       nodes_local(2,2) = n2
       nodes_local(3,3) = n3
       nodes_local(4,4) = n4
       nodes_local(2,1) = n12
       nodes_local(3,1) = n13
       nodes_local(4,1) = n14
       nodes_local(3,2) = n23
       nodes_local(4,2) = n24
       nodes_local(4,3) = n34
       ! build first four tetrahedrons
       this%topol(1,ntetref(1)) = n1  ; this%topol(1,ntetref(2)) = n2  ; this%topol(1,ntetref(3)) = n3  ; this%topol(1,ntetref(4)) = n4  ; 
       this%topol(2,ntetref(1)) = n12 ; this%topol(2,ntetref(2)) = n12 ; this%topol(2,ntetref(3)) = n13 ; this%topol(2,ntetref(4)) = n14 ; 
       this%topol(3,ntetref(1)) = n13 ; this%topol(3,ntetref(2)) = n23 ; this%topol(3,ntetref(3)) = n23 ; this%topol(3,ntetref(4)) = n24 ; 
       this%topol(4,ntetref(1)) = n14 ; this%topol(4,ntetref(2)) = n24 ; this%topol(4,ntetref(3)) = n34 ; this%topol(4,ntetref(4)) = n34 ;
       ! evaluate Euclidean norm between opposite middle edge nodes
       ! and select the case with minimum distance
       difference1_x=this%coord(1,n12)-this%coord(1,n34)
       difference1_y=this%coord(2,n12)-this%coord(2,n34)
       difference1_z=this%coord(3,n12)-this%coord(3,n34)
       difference2_x=this%coord(1,n13)-this%coord(1,n24)
       difference2_y=this%coord(2,n13)-this%coord(2,n24)
       difference2_z=this%coord(3,n13)-this%coord(3,n24)
       difference3_x=this%coord(1,n14)-this%coord(1,n23)
       difference3_y=this%coord(2,n14)-this%coord(2,n23)
       difference3_z=this%coord(3,n14)-this%coord(3,n23)
       !
       norm1=sqrt((difference1_x*difference1_x)+(difference1_y*difference1_y)+(difference1_z*difference1_z))
       norm2=sqrt((difference2_x*difference2_x)+(difference2_y*difference2_y)+(difference2_z*difference2_z))
       norm3=sqrt((difference3_x*difference3_x)+(difference3_y*difference3_y)+(difference3_z*difference3_z))
       ! build other four tetrahedrons, whose construction
       ! depends on the norm evaluated above
       if ((norm1 .le. norm2) .and. (norm1 .le. norm3)) then
          ! case 1, minimum distance between n12 and n34
          ! the first two nodes are n12 and n34
          this%topol(1,ntetref(5):ntetref(8)) = n12
          this%topol(2,ntetref(5):ntetref(8)) = n34
          this%topol(3,ntetref(5):ntetref(6)) = n14
          this%topol(3,ntetref(7):ntetref(8)) = n23
          this%topol(4,ntetref(5))            = n24
          this%topol(4,ntetref(6))            = n13
          this%topol(4,ntetref(7))            = n24
          this%topol(4,ntetref(8))            = n13
       elseif ((norm2 .le. norm1) .and. (norm2 .le. norm3)) then
          ! case 2, minimum distance between n13 and n24
          ! the first two nodes are n13 and n24
          this%topol(1,ntetref(5):ntetref(8)) = n13
          this%topol(2,ntetref(5):ntetref(8)) = n24
          this%topol(3,ntetref(5):ntetref(6)) = n34
          this%topol(3,ntetref(7):ntetref(8)) = n12
          this%topol(4,ntetref(5))            = n14
          this%topol(4,ntetref(6))            = n23
          this%topol(4,ntetref(7))            = n14
          this%topol(4,ntetref(8))            = n23
       elseif ((norm3 .le. norm1) .and. (norm3 .le. norm2)) then
          ! case 3, minimum distance between n14 and n23
          ! the first two nodes are n14 and n23
          this%topol(1,ntetref(5):ntetref(8)) = n14
          this%topol(2,ntetref(5):ntetref(8)) = n23
          this%topol(3,ntetref(5):ntetref(6)) = n12
          this%topol(3,ntetref(7):ntetref(8)) = n34
          this%topol(4,ntetref(5))            = n24
          this%topol(4,ntetref(6))            = n13
          this%topol(4,ntetref(7))            = n24
          this%topol(4,ntetref(8))            = n13
       end if
       
       do iloc=1,8
          this%topol(5,ntetref(iloc)) = izone
          this%cell_parent(ntetref(iloc)) = icell
          ! sort the nodes of refined mesh in increasing order
          call isort(4,this%topol(1:4,ntetref(iloc)))
       end do

       ! update column numbers of the refined tetrahedrons
       ! in the array this%topol 
       ntetref(:) = ntetref(:) + 8
       
    end do
    
    call local_mesh%kill(lun_err)
   
  end subroutine build_refined_topol_coord

  !>-------------------------------------------------------------
  !> Static destructor.
  !> (procedure public for type mesh)
  !> deallocate all arrays for a var of type mesh
  !>
  !> usage:
  !>     call 'var'%kill(lun)
  !>
  !> where:
  !> \param[in] lun -> integer. I/O unit for error message output
  !>
  !<-----------------------------------------------------------
  subroutine kill_mesh(this, lun)
    implicit none
    ! vars
    class(mesh3d) :: this
    integer, intent(in) :: lun
    ! local vars
    integer :: res
    logical :: rc
    integer :: inode, icell
    
    if (allocated(this%topol)) then
       deallocate(this%topol,stat=res)
       if (res.ne.0) rc=IOerr(lun, err_dealloc, 'kill_mesh', &
            'mesh members topol')
    end if

    if (allocated(this%coord)) then
       deallocate(this%coord,stat=res)
       if (res.ne.0) rc=IOerr(lun, err_dealloc, 'kill_mesh', &
            'mesh members coord')
    end if

    if (allocated(this%size_cell)) then
       deallocate(this%size_cell,stat=res)
       if (res.ne.0) rc=IOerr(lun, err_dealloc, 'kill_mesh', &
            'mesh members size_cell')
    end if

    if (allocated(this%size_cell)) then
      deallocate(this%size_node,stat=res)
       if (res.ne.0) rc=IOerr(lun, err_dealloc, 'kill_mesh', &
            'mesh members size_node')
    end if
    
       
    if (allocated(this%bar_cell)) then
       deallocate(this%bar_cell,stat=res)
       if (res.ne.0) rc=IOerr(lun, err_dealloc, 'kill_mesh', &
            'mesh members bar_cell')
    end if

    ! the followings are done just in case this%cnc_built greater than 0
    
    if (allocated(this%surface_tetra)) then
       deallocate(this%surface_tetra,stat=res)
       if (res.ne.0) rc=IOerr(lun, err_dealloc, 'kill_mesh', &
            'mesh members surface_tetra')
    end if

    if (allocated(this%normal)) then
       deallocate(this%normal,stat=res)
       if (res.ne.0) rc=IOerr(lun, err_dealloc, 'kill_mesh', &
            'mesh members normal')
    end if

    if (allocated(this%surface_face)) then
       deallocate(this%surface_face,stat=res)
       if (res.ne.0) rc=IOerr(lun, err_dealloc, 'kill_mesh', &
            'mesh members surface_face')
    end if

    if (allocated(this%bar_edge)) then
       deallocate(this%bar_edge,stat=res)
       if (res.ne.0) rc=IOerr(lun, err_dealloc, 'kill_mesh', &
            'mesh members bar_edge')
    end if

    if (allocated(this%face_plist)) then
       deallocate(this%face_plist,stat=res)
       if (res.ne.0) rc=IOerr(lun, err_dealloc, 'kill_mesh', &
            'mesh members face_plist')
    end if

    if (allocated(this%face_cnc)) then
       deallocate(this%face_cnc,stat=res)
       if (res.ne.0) rc=IOerr(lun, err_dealloc, 'kill_mesh', &
            'mesh members face_cnc')
    end if

    if (allocated(this%face_iside)) then
       deallocate(this%face_iside,stat=res)
       if (res.ne.0) rc=IOerr(lun, err_dealloc, 'kill_mesh', &
            'mesh members face_iside')
    end if

    if (allocated(this%neigh_face)) then
       deallocate(this%neigh_face,stat=res)
       if (res.ne.0) rc=IOerr(lun, err_dealloc, 'kill_mesh', &
            'mesh members neigh_face')
    end if

    if (allocated(this%edge_cnc)) then
       deallocate(this%edge_cnc,stat=res)
       if (res.ne.0) rc=IOerr(lun, err_dealloc, 'kill_mesh', &
            'mesh members edge_cnc')
    end if

    if (allocated(this%iside_edge)) then
       deallocate(this%iside_edge,stat=res)
       if (res.ne.0) rc=IOerr(lun, err_dealloc, 'kill_mesh', &
            'mesh members iside_edge')
    end if

    this%nnode = 0
    this%ncell = 0

    ! this is done just in case this%grid_level greater than 0
    if (allocated(this%cell_parent)) then
       deallocate(this%cell_parent,stat=res)
       if (res.ne.0) rc=IOerr(lun, err_dealloc, 'kill_mesh', &
            'mesh members cell_parent')
    end if

    if (allocated(this%node_parent)) then
       deallocate(this%node_parent,stat=res)
       if (res.ne.0) rc=IOerr(lun, err_dealloc, 'kill_mesh', &
            'mesh members node_parent')
    end if
       
    this%grid_level=0

    ! this is done if flag_reorder is present and is equal to 1
    if (allocated(this%perm_nodes)) then
       deallocate(this%perm_nodes,stat=res)
       if (res.ne.0) rc=IOerr(lun, err_dealloc, 'kill_mesh', &
            'mesh members perm_nodes')
    end if
    
  end subroutine kill_mesh
  !>-------------------------------------------------------------
  !> Info procedure.
  !> (public procedure for type mesh)
  !> Prints content of a variable of type mesh
  !> 
  !> usage:
  !>     call 'var'%info(lun)
  !>
  !> where:
  !> \param[in] lun -> integer. I/O unit for error message output
  !>
  !<-------------------------------------------------------------
  subroutine info_mesh(this, lun)
    use Globals
    implicit none
    ! vars
    class(mesh3d), intent(in) :: this
    integer , intent(in) :: lun
    ! local vars
    integer :: j, k
    
    write(lun,*) ' '
    write(lun,*) ' Info: grid structure definition:'

    write(lun,*) 'number of nodes        ',this%nnode
    write(lun,*) 'number of tetrahedrons ',this%ncell
    write(lun,*) 'number of total faces  ',this%nface
    write(lun,*) 'number of bnd faces    ',this%nface_bc
    write(lun,*) 'number of int faces    ',this%nface-this%nface_bc
    write(lun,*) 'number of edges        ',this%nedge
    write(lun,*) 'number of bnd nodes   ',this%nnode_bc

    write(lun,*) ' '
    write(lun,*) 'node coordinates: (type mesh member coord)'
    write(lun,*) '    node   x           y           z'
    do j=1,this%nnode
       write(lun,'(i8,2x,3(1pe12.4))') j,(this%coord(k,j),k=1,3)
    end do

    write(lun,*) ' '
    write(lun,*) 'tetrahedrons: (type mesh memeber Tetrahedron)'
    write(lun,*) '    tetra    node 1    node 2    node 3    node 4   material'
    do j=1,this%ncell
       write(lun,'(i8,2x,i8,2x,i8,2x,i8,2x,i8,2x,i5)') j,(this%topol(k,j),k=1,5)
    end do

    write(lun,*) ' '    
    write(lun,*) 'element volume: (type mesh member volume)'
    do j=1,this%ncell
       write(lun,'(i8,2x,(1pe12.4))') j,this%size_cell(j)
    end do

    write(lun,*) ' '    
    write(lun,*) 'nodal Volume: (type mesh member size_node)'
    do j=1,this%nnode
       write(lun,'(i8,2x,(1pe12.4))') j,this%size_node(j)
    end do


    if (this%cnc_built.ge.1) then
       write(lun,*) ' '    
       write(lun,*) 'Boundary faces: (type mesh member face_iside)'
       write(lun,*) '   face    node 1    node 2    node 3   '
       do j=1,this%nface_bc
       write(lun,'(i8,2x,i8,2x,i8,2x,i8)') j,(this%face_iside(k,j),k=1,3)
       end do

       write(lun,*) ' '    
       write(lun,*) 'Internal faces: (type mesh member face_iside)'
       write(lun,*) '   face    node 1    node 2   node 3'
       do j=this%nface_bc+1,this%nface
       write(lun,'(i8,2x,i8,2x,i8,2x,i8)') j,(this%face_iside(k,j),k=1,3)
       end do
       
       write(lun,*) ' '    
       write(lun,*) 'neighbors: (type mesh member neigh_face)'
       write(lun,*) '    tetra     tetra 1     tetra 2     tetra 3     tetra 4 '
       do j=1,this%ncell
          write(lun,'(i8,2x,i8,2x,i8,2x,i8,2x,i8)') j,(this%neigh_face(k,j),k=1,4)
       end do
       
       write(lun,*) ' '    
       write(lun,*) 'face numbering: (type mesh member face_cnc)'
       write(lun,*) '    tetra    face 1    face 2    face 3   face 4 '
       do j=1,this%ncell
          write(lun,'(i8,2x,i8,2x,i8,2x,i8,2x,i8)') j,(this%face_cnc(k,j),k=1,4)
       end do
       
       write(lun,*) ' '    
       write(lun,*) 'west-east elements: (type mesh member face_plist)'
       write(lun,*) '   face     tetra 1     tetra 2 '
       do j=1,this%nface
          write(lun,'(i8,2x,i8,2x,i8)') j,(this%face_plist(k,j),k=1,2)
       end do

       write(lun,*) ' '    
       write(lun,*) 'element surface: (type mesh member surface)'
       do j=1,this%ncell
          write(lun,'(i8,2x,(1pe12.4))') j,this%surface_tetra(j)
       end do

       write(lun,*) ' '    
       write(lun,*) 'surface normal : (type mesh member normal) (non ancora analizzato il verso nelle facce sul bordo)'
       write(lun,*) '   face    s1n       s2n       s3n '
       do j=1,this%nface
          write(lun,'(i8,2x,3(1pe12.4))') j,(this%normal(k,j),k=1,3)
       end do

       
       write(lun,*) ' '    
       write(lun,*) 'coordinates of cell baricenters: (type mesh members bar_cell)'
       write(lun,*) '   elem     s1c     s2c     s3c     '
       do j=1,this%ncell
          write(lun,'(i8,2x,3(1pe12.4))') j,(this%bar_cell(k,j),k=1,3)
       end do

       write(lun,*) ' '    
       write(lun,*) 'Edges: (type mesh member iside_edge)'
       write(lun,*) '   edge    node 1    node 2'
       do j=1,this%nedge
       write(lun,'(i8,2x,i8,2x,i8)') j,(this%iside_edge(k,j),k=1,2)
       end do

       write(lun,*) ' '    
       write(lun,*) 'edge numbering: (type mesh member edge_cnc)'
       write(lun,*) '    tetra    edge 1    edge 2    edge 3   edge 4    edge 5    edge 6 '
       do j=1,this%ncell
          write(lun,'(i8,2x,i8,2x,i8,2x,i8,2x,i8,2x,i8,2x,i8)') j,(this%edge_cnc(k,j),k=1,6)
       end do
       
    end if
    
    if (this%grid_level.gt.0) then
       write(lun,*) ' '    
       write(lun,*) 'parent information'
       write(lun,*) ' '    
       write(lun,*) ' node_parent'    
       write(lun,*) '   node     par 1    par 2'
       do j=1,this%nnode
          write(lun,'(i8,2x,i8,2x,i8)') j,(this%node_parent(k,j),k=1,2)
       end do
       write(lun,*) ' '    
       write(lun,*) ' cell_parent'    
       write(lun,*) '   tria     parent'
       do j=1,this%ncell
          write(lun,'(i8,2x,i8)') j,this%cell_parent(j)
       end do
    end if
  end subroutine info_mesh

  !>--------------------------------------------------------------
  !> Procedure write the mesh in the format
  !> nnode
  !> nncell
  !> coord1,coord2, coord3
  !> node1,node2,node3,node4
  !>-----------------------------------------------------------
  subroutine write_mesh(this, lun_err, file2write)
    use Globals
    implicit none
    ! vars
    class(mesh3d), intent(in) :: this
    integer,     intent(in) :: lun_err
    type(file),  intent(in) :: file2write
    ! local vars
    logical :: rc
    integer :: res,lun_write
    integer :: inode, icell,j
    
    lun_write=file2write%lun
    write(lun_write,*,iostat=res) this%nnode, '   #  of nodes '
    if(res .ne. 0) rc = IOerr(lun_err, err_out, 'write mesh', &
            'type mesh member nnode',res)

    write(lun_write,*,iostat=res) this%ncell, this%nnodeincell, ' 3d ! # of cells,  # of node in cell,  mesh type'
    if(res .ne. 0) rc = IOerr(lun_err, err_out, 'write mesh', &
         'type mesh member ncell',res)

    do inode = 1,this%nnode
       write(lun_write,*,iostat=res) this%coord(1:3,inode)
       if(res .ne. 0) rc = IOerr(lun_err, err_out, 'write mesh', &
            'type mesh meber coord',res)
    end do
    do icell = 1, this%ncell
       write(lun_write,*,iostat=res) ( this%topol(j,icell), j=1,5 )
       if(res .ne. 0) rc = IOerr(lun_err, err_out, 'write mesh', &
            'type mesh member tetra',res)
    end do
  end subroutine write_mesh

  !>-------------------------------------------------------------
  !> calculate the projection along the edge outward unit normal
  !> (public procedure for type mesh)
  !>
  !>
  !> usage:
  !>     call 'var'%proj(vec, proj)
  !>
  !> where:
  !> \param[in] vec -> real(double), dimension (2,nedge)
  !>                   input vector to be projected along the edge normal
  !> \param[out] proj-> real(double), dimension (nedge)
  !>                    array of projected velocity along the normal to the edge
  !>                    It is a scalar for each edge as it is assumed to
  !>                    be directed along the normal.
  !>
  !<-------------------------------------------------------------
  subroutine normal_projection(this, vec, proj)
    use Globals
    implicit none

    class(mesh3d), intent(in) :: this
    real(kind=double), intent(in) :: vec(3,this%nface_bc)
    real(kind=double), intent(out) :: proj(this%nface_bc)

    ! local vars
    integer :: iface
    real(kind=double) :: locvec(3),locnorm(3)
    ! external functions (BLAS lvl 1)
    real(kind=double) :: ddot

    do iface=1,this%nface_bc
       locvec = vec(:,iface)
       locnorm = this%normal(:,iface)
       proj(iface)=ddot(3,locvec,1,locnorm,1)
    end do
  end subroutine normal_projection
  !>-------------------------------------------------------------
  !> Build level 1 mesh connectivity data starting from 
  !> previous mesh%topol mesh%coord
  !> (private procedure for type mesh, used in init)
  !> 
  !> usage:
  !>     call 'var'%connection_faces(lun)
  !>
  !> where:
  !> \param[in] lun -> integer. I/O unit for error message output
  !>
  !> part of the static constructor
  !> at the end, we permutate all the face structure and we
  !> renumber the faces so that the first nface_bc faces are at
  !> the boundary
  !<-------------------------------------------------------------
  !! we work with local copies, and do proper assignments
  !! at the beginning and at the end
  subroutine connection_faces(this, lun)
    use Globals
    implicit none
    !  External variables 
    class(mesh3d), intent(inout) :: this
    integer,     intent(in) :: lun
    ! Local variables
    ! work arrays
    integer, allocatable :: topol(:,:)
    integer, allocatable :: faces(:,:), point(:)
    integer, allocatable :: neigh(:,:),side_cnc(:,:)
    integer, allocatable :: iside(:,:),plist(:,:)
    integer, allocatable :: perm_faces(:),perminv_faces(:)
    integer, allocatable :: mark_node(:),loc_node_bc(:)
    !  scalars and static arrays
    integer :: perm(3,4)
    data perm/1,2,3, 1,2,4, 1,3,4, 2,3,4/
    integer :: nface,nface_bc,ncell,nnode
    integer :: res
    integer :: icell,jtetra,ktetra,jj,kk,iloc,inode
    integer :: iface,jface,kface
    integer :: i,j,m,indx,isgn, itemp
    integer :: temp(4)
    logical :: rc


    ! Assign local scalar variabls
    ncell = this%ncell
    nnode = this%nnode
    !
    ! Allocate local and work arrays
    ! 
    allocate(topol(4,ncell),faces(5,4*ncell),point(4*ncell),side_cnc(4,ncell),plist(2,4*ncell),iside(3,4*ncell),neigh(4,ncell),stat=res)
    if(res .ne. 0) rc = IOerr(lun, err_alloc, 'connection_faces', &
         'work and local arrays topol, faces, point, side_cnc, plist, iside, neigh',res)
    !
    ! generate faces vector
    !
    topol=this%topol(1:4,:)
    do icell = 1,ncell
       !
       ! sort the element nodes in increasing order
       !
       call isort(4,topol(1:4,icell))
       do jtetra = 1,4
          kk = 4*(icell-1) + jtetra
          faces(1,kk) = topol(perm(1,jtetra),icell)
          faces(2,kk) = topol(perm(2,jtetra),icell)
          faces(3,kk) = topol(perm(3,jtetra),icell)
          faces(4,kk) = icell
          faces(5,kk) = 0
       enddo
    enddo
    !
    ! sort in lexicographic order the faces vector
    !
    !call heapsort_faces(faces,4*ncell)
    
    !  Initialize.
    i = 0
    indx = 0
    isgn = 0
    j = 0
    do 
       call global_heapsort(4*ncell, indx, i,j,isgn)
       if (indx .gt. 0 ) then
          ! SWAP ELEMENT 

          ! swap first 4 element of faces
          temp(1:4)    = faces(1:4,i)
          faces(1:4,i) = faces(1:4,j)
          faces(1:4,j) = temp

       else if ( indx .lt. 0) then
          ! COMPARE
          isgn = 0
          if ( lexicographic_order( 3,faces(1:3,i),faces(1:3,j) ) ) then
             isgn = -1
          else
             isgn = 1
          end if
       else if ( indx .eq. 0 ) then
          exit
       end if
    end do

    !
    ! eliminate double faces and build side_cnc
    !
    nface = 1
    do icell = 2,4*ncell
       !
       ! Check the faces, if coincident remove the newer
       ! and update faces(5,nface) so that faces(4,nface) < faces(5,nface)
       if ((faces(1,icell) .eq. faces(1,nface)) .and.                      &
            &       (faces(2,icell) .eq. faces(2,nface)) .and.                      &
            &       (faces(3,icell) .eq. faces(3,nface))) then
          if (faces(4,icell) .gt. faces(4,nface)) then
             faces(5,nface) = faces(4,icell)
          else 
             faces(5,nface) = faces(4,nface)
             faces(4,nface) = faces(4,icell)
          end if
       else
          nface = nface + 1
          do jtetra = 1,5
             faces(jtetra,nface) = faces(jtetra,icell)
          end do
       end if
    end do
    !
    !  build the neigh vector using plist as a work array
    !
    do icell = 1,nface
       point(icell) = 1
    end do

      ! set side_cnc and neigh to 0
      side_cnc = 0
      neigh = 0

      do icell = 1,nface
         iside(1,icell) = faces(1,icell)
         iside(2,icell) = faces(2,icell)
         iside(3,icell) = faces(3,icell)
         plist(1,icell) = faces(4,icell)
         plist(2,icell) = faces(5,icell)
         jj = faces(4,icell)
         kk = faces(5,icell)
         if (kk .ne. 0) then

!           face inside the domain
            neigh(point(jj),jj) = kk
            neigh(point(kk),kk) = jj
            point(jj) = point(jj) + 1
            point(kk) = point(kk) + 1            
         end if
      end do

      do icell = 1,nface
         jtetra = faces(4,icell)
         ktetra = faces(5,icell)
         kk = 1
         do while (kk .le. 4)
            if (side_cnc(kk,jtetra) .eq. 0) then
               side_cnc(kk,jtetra) = icell
               kk = 5
            end if
            kk = kk + 1
         end do
         if (ktetra .gt. 0) then
            kk = 1
            do while (kk .le. 4)
               if (side_cnc(kk,ktetra) .eq. 0) then
                  side_cnc(kk,ktetra) = icell
                  kk = 5
               end if
               kk = kk + 1
            end do
         end if
      end do

      allocate(perm_faces(nface),perminv_faces(nface),stat=res)
      if(res .ne. 0) rc = IOerr(lun, err_alloc, 'connection_faces', &
           ' temporary array perm_faces, perminv_faces ', res)
     
      kface=0
      do iface=1,nface
         if(plist(1,iface)*plist(2,iface).eq.0) then
            ! for the faces at the boundary
            kface=kface+1
            perm_faces(kface)=iface
         end if
      end do
      nface_bc=kface
      do iface=1,nface
         if(plist(1,iface)*plist(2,iface).ne.0) then
            ! for the internal faces
            kface=kface+1
            perm_faces(kface)=iface
         end if
      end do
      do iface=1,nface
         perminv_faces(perm_faces(iface))=iface
      end do
      !
      ! Allocate global arrays of the structure "mesh"
      ! if not yet allocated
      !
      if ( .not. allocated(this%neigh_face) ) then
         allocate(this%neigh_face(4,ncell),stat=res)
         if(res .ne. 0) rc = IOerr(lun, err_alloc, 'connection_faces', &
              ' type mesh member array neigh_face',res)
      end if
      if ( .not. allocated(this%side_cnc) ) then
         allocate(this%side_cnc(4,ncell),stat=res)
         if(res .ne. 0) rc = IOerr(lun, err_alloc, 'connection_faces', &
              ' type mesh member array side_cnc',res)
      end if

      if ( .not. allocated(this%iside_face) ) then
         allocate(this%iside_face(3,nface),stat=res)
         if(res .ne. 0) rc = IOerr(lun, err_alloc, 'connection_faces', &
              ' type mesh member array iside_face',res)
      end if
      if ( .not. allocated(this%plist) ) then
         allocate(this%plist(2,nface),stat=res)
         if(res .ne. 0) rc = IOerr(lun, err_alloc, 'connection_faces', &
              ' type mesh member array plist',res)
      end if

      this%nface = nface
      this%nface_bc = nface_bc

      do iface=1,nface
         this%iside_face(:,iface)=iside(:,perm_faces(iface))
         this%plist(:,iface)=plist(:,perm_faces(iface))
      end do
      do icell=1,ncell
         do iloc=1,4
            this%side_cnc(iloc,icell)=perminv_faces(side_cnc(iloc,icell))
         end do
      end do
      
      this%neigh_face=neigh

      ! boundary_nodes
      allocate(&
           loc_node_bc(nnode),&
           mark_node(nnode),&
           stat=res)
      if(res .ne. 0) rc = IOerr(lun, err_dealloc, 'connection_faces', &
           ' local arrays mark_node',res)
      mark_node = 0
      do iface=1,this%nface_bc
         mark_node(this%iside_face(1:3,iface)) = 1
      end do
      m=0
      do inode=1,nnode
         if( mark_node(inode) .eq. 1 ) then
            m=m+1
            loc_node_bc(m) = inode
         end if
      end do
      allocate(&
           this%node_bc(m),&
           stat=res)
      if(res .ne. 0) rc = IOerr(lun, err_dealloc, 'connection_faces', &
           ' array node_bc',res)
      this%nnode_bc = m
      this%node_bc = loc_node_bc(1:m)
      
    ! deallocate local and work arrays
    deallocate(faces,point,topol,side_cnc,plist,iside,neigh,perm_faces,perminv_faces,mark_node,loc_node_bc,stat=res)
    if(res .ne. 0) rc = IOerr(lun, err_dealloc, 'connection_faces', &
    'work arrays faces, point, topol, side_cnc, plist, iside, neigh, perm_faces, perminv_faces, perm, mark_node,loc_node_bc',res)

  contains

!!$    subroutine heapsort_faces(face,n_aux)
!!$  !------------------------------------------------------------
!!$  ! sort array face of integers in ascending numerical order
!!$  !-------------------------------------------------------------
!!$      use Globals
!!$      implicit none
!!$      integer  node,i_aux,j_aux,k_aux,ik,jk
!!$      integer  n_aux 
!!$      integer  face(5,n_aux)
!!$
!!$      do node = 2,n_aux
!!$         i_aux = node
!!$         j_aux = i_aux/2
!!$         do while((i_aux.ne.1).and.(conf_order(3,face(1:3,j_aux),face(1:3,i_aux))))      
!!$            call swap_faces(face(1,j_aux),face(1,i_aux))
!!$            i_aux = j_aux
!!$            j_aux = i_aux/2
!!$         end do
!!$      end do
!!$      do i_aux = n_aux,2,-1
!!$         call swap_faces(face(1,i_aux),face(1,1))
!!$         k_aux = i_aux-1
!!$         ik = 1
!!$         jk = 2
!!$
!!$         if ((k_aux.ge.3).and.(conf_order(3,face(1:3,2),face(1:3,3))))  jk = 3
!!$         do while ((jk.le.k_aux).and.(conf_order(3,face(1:3,ik),face(1:3,jk))))
!!$            call swap_faces(face(1,jk),face(1,ik))
!!$            ik = jk
!!$            jk = ik*2
!!$            if ((jk+1 .le. k_aux) .and. (conf_order(3,face(1:3,jk),face(1:3,jk+1))))&
!!$     &         jk = jk+1  
!!$
!!$         end do
!!$      end do
!!$      return
!!$    end subroutine heapsort_faces
    
    subroutine swap_faces(a_aux,b_aux)
      implicit none
      integer  a_aux(4),b_aux(4)
      integer  i_aux,temp
      do i_aux=1,4
         temp = a_aux(i_aux)
         a_aux(i_aux) = b_aux(i_aux)
         b_aux(i_aux) = temp
      end do
      return
    end subroutine swap_faces

  end subroutine connection_faces

  !>-------------------------------------------------------------
  !> Build level 1 mesh edges connectivity data starting from 
  !> previous mesh%topol mesh%coord
  !> (private procedure for type mesh, used in init)
  !> 
  !> usage:
  !>     call 'var'%connection_edges(lun)
  !>
  !> where:
  !> \param[in] lun -> integer. I/O unit for error message output
  !>
  !<-------------------------------------------------------------
  !! we work with local copies, and do proper assignments
  !! at the beginning and at the end
  subroutine connection_edges(this, lun)
    use Globals
    implicit none
    !  External variables 
    class(mesh3d), intent(inout) :: this
    integer,     intent(in   ) :: lun
    ! Local variables
    ! work arrays
    integer, allocatable :: topol(:,:)
    integer, allocatable :: edges(:,:), point(:)
    integer, allocatable :: neigh(:,:),edge_cnc(:,:)
    integer, allocatable :: iside(:,:),plist(:,:)
    integer, allocatable :: perm_faces(:),perminv_faces(:)
    !  scalars and static arrays
    integer :: perm(2,6)
    data perm/1,2, 1,3, 1,4, 2,3, 2,4, 3,4/
    integer :: ncell, nedge, count_edge
    integer :: res
    integer :: icell,jtetra,ktetra,iedge,jedge,jj,kk,iloc
    integer :: i,j, indx,isgn, itemp
    integer :: temp(3)
    logical :: rc

    ! Assign local scalar variables
    ncell=this%ncell
    nedge=6*ncell
    !
    ! Allocate local and work arrays
    !
    allocate(topol(5,ncell),edges(4,nedge),point(ncell),edge_cnc(6,ncell),iside(2,6*ncell),stat=res)
    if(res .ne. 0) rc = IOerr(lun, err_alloc, 'connection_edges', &
         'work and local arrays topol, edges, point, edge_cnc, iside',res)
    !
    ! Generate edges vector
    !
    topol=this%topol

    do icell = 1,ncell
       ! sort the element nodes in increasing order
       call isort(4,topol(1:4,icell))
       do jtetra = 1,6
          kk = 6*(icell-1) + jtetra
          edges(1,kk) = topol(perm(1,jtetra),icell)
          edges(2,kk) = topol(perm(2,jtetra),icell)
          edges(3,kk) = icell
          edges(4,kk) = 0
       enddo
    enddo
    !
    ! sort in lexicographic order the edges vector
    !
    
    ! call heapsort_edges(edges,nedge)
    !  Initialize.
    i = 0
    indx = 0
    isgn = 0
    j = 0
    do 
       call global_heapsort(nedge, indx, i,j,isgn)
       if (indx .gt. 0 ) then
          ! SWAP ELEMENT 

          ! swap first 3 element of edges
          temp(1:3)    = edges(1:3,i)
          edges(1:3,i) = edges(1:3,j)
          edges(1:3,j) = temp

       else if ( indx .lt. 0) then
          ! COMPARE
          isgn = 0
          if ( lexicographic_order( 2,edges(1:2,i),edges(1:2,j) ) ) then
             isgn = -1
          else
             isgn = 1
          end if
       else if ( indx .eq. 0 ) then
          exit
       end if
    end do
    
    


!!$    do j=1,nedge
!!$       write(*,'(4I2)') (edges(i,j),i=1,4)
!!$    end do
    ! NOTA CHE NON ORDINA TROPPO PERFETTAMENTE, IN PARTICOLARE LA TERZA COLONNA NE SBAGLIA UNO
    !
    ! update number of node
    !
    count_edge = 1
    edges(4,1) = count_edge
    do iedge=2,nedge
       if ((edges(1,iedge) .eq. edges(1,iedge-1)) .and. &
            (edges(2,iedge) .eq. edges(2,iedge-1))) then
          edges(4,iedge) = count_edge
       else
          count_edge = count_edge + 1
          edges(4,iedge) = count_edge
       end if
    end do
    !
    ! build iside
    !
    iside(1:2,1) = edges(1:2,1)
    count_edge = 1
    do iedge=2,nedge
       if (edges(4,iedge) .ne. count_edge) then
          count_edge = count_edge + 1
          iside(:,count_edge) = edges(1:2,iedge)
       end if
    end do
    this%nedge = count_edge
    !
    ! build edge_cnc
    ! initialize edge_cnc to zero
    !
    edge_cnc = 0
    do iedge=1,nedge
       jtetra = edges(3,iedge)
       kk = 1
       do while (kk .le. 6)
          if (edge_cnc(kk,jtetra) .eq. 0) then
             edge_cnc(kk,jtetra) = edges(4,iedge)
             kk = 7
          end if
          kk = kk + 1
       end do
    end do
    !
    ! Allocate global arrays of the structure "mesh"
    ! if not yet allocated
    !
    if ( .not. allocated(this%edge_cnc) ) then
       allocate(this%edge_cnc(6,this%ncell),stat=res)
       if(res .ne. 0) rc = IOerr(lun, err_alloc, 'connection_edges', &
            ' type mesh member array edge_cnc',res)
    end if
    if ( .not. allocated(this%iside_edge) ) then
       allocate(this%iside_edge(2,this%nedge),stat=res)
       if(res .ne. 0) rc = IOerr(lun, err_alloc, 'connection_edges', &
            ' type mesh member array iside_edge',res)
    end if

    this%edge_cnc = edge_cnc
    this%iside_edge = iside

    
    ! deallocate local and work arrays
    deallocate(edges,point,topol,edge_cnc,iside,stat=res)
    if(res .ne. 0) rc = IOerr(lun, err_dealloc, 'connection_edges', &
    'work arrays edges, point, topol, edge_cnc, iside',res)
    
    contains

!!$      subroutine heapsort_edges(face,n_aux)
!!$        !>---------------------------------------------------------
!!$        !> sort array edge of integers in ascending numerical order
!!$        !>---------------------------------------------------------
!!$        use Globals
!!$        implicit none
!!$        integer  node,i_aux,j_aux,k_aux,ik,jk
!!$        integer  n_aux 
!!$        integer  face(4,n_aux)      
!!$        do node = 2,n_aux
!!$           i_aux = node
!!$           j_aux = i_aux/2
!!$           do while((i_aux.ne.1).and.(conf_order(2,face(1:2,j_aux),face(1:2,i_aux))))
!!$              call swap_edges(face(1:3,j_aux),face(1:3,i_aux))
!!$            i_aux = j_aux
!!$            j_aux = i_aux/2
!!$         end do
!!$      end do
!!$      do i_aux = n_aux,2,-1
!!$         call swap_edges(face(1:3,i_aux),face(1:3,1))
!!$         k_aux = i_aux-1
!!$         ik = 1
!!$         jk = 2
!!$         if ((k_aux.ge.3).and.(conf_order(2,face(1:2,2),face(1:2,3))))  jk = 3
!!$         do while ((jk.le.k_aux).and.(conf_order(2,face(1:2,ik),face(1:2,jk))))
!!$            call swap_edges(face(1:3,jk),face(1:3,ik))
!!$            ik = jk
!!$            jk = ik*2 
!!$            if ((jk+1 .le. k_aux) .and. (conf_order(2,face(1:2,jk),face(1:2,jk+1))))&
!!$     &         jk = jk+1  
!!$            
!!$         end do
!!$      end do
!!$    end subroutine heapsort_edges
    
    subroutine swap_edges(a_aux,b_aux)
      implicit none
      integer  a_aux(3),b_aux(3)
      integer  i_aux,temp
      do i_aux=1,3
         temp = a_aux(i_aux)
         a_aux(i_aux) = b_aux(i_aux)
         b_aux(i_aux) = temp
      end do
    end subroutine swap_edges
      
  end subroutine connection_edges

  
    !>-------------------------------------------------------------
    !> calculate cell volume, edge, lenght, normals, etc
    !> completing the instantiation of a var of type mesh
    !> (private procedure for type mesh, used in init)
    !> 
    !> usage:
    !>     call 'var'%meshProp(lun)
    !>
    !> where:
    !> \param[in] lun -> integer. I/O unit for error message output
    !>
    !> part of the static constructor
    !<-------------------------------------------------------------
    subroutine meshProp(this, lun)
      use Globals
      implicit none
      ! vars
      class(mesh3d), intent(inout) :: this
      integer, intent(in)        :: lun    
      ! local variables
      integer           :: n1,n2,n3,n4    
      real(kind=double) :: s1n1,s2n1,s3n1,s1n2,s2n2,s3n2,s1n3,s2n3,s3n3, s1n4,s2n4,s3n4
      real(kind=double) :: r1,r2,r3,d1,d2,d3
      integer           :: icell,iface,iloc,iedg,inode
      integer           :: res
      logical           :: rc

      if (.not. allocated(this%size_cell)) then
         allocate(this%size_cell(this%ncell),stat=res)
         if(res .ne. 0) rc = IOerr(lun, err_alloc, 'meshProp', &
              '  type mesh member volume (array)',res)
      end if
      if (.not. allocated(this%size_node)) then
         allocate(this%size_node(this%nnode),stat=res)
         if(res .ne. 0) rc = IOerr(lun, err_alloc, 'meshProp', &
              '  type mesh member size_node',res)
      end if

      if(.not. allocated(this%bar_cell)) then
         allocate(this%bar_cell(3,this%ncell),stat=res)
         if(res .ne. 0) rc = IOerr(lun, err_alloc, 'meshProp', &
              '  type mesh member bar_cell (array)',res)
      end if
      if (this%cnc_built .ge. 1) then
          if (.not. allocated(this%surface_face)) then
            allocate(this%surface_face(this%nface),stat=res)
            if(res .ne. 0) rc = IOerr(lun, err_alloc, 'meshProp', &
                 '  type mesh member surface_face (array)',res)
         end if
         if (.not. allocated(this%surface_tetra)) then
            allocate(this%surface_tetra(this%ncell),stat=res)
            if(res .ne. 0) rc = IOerr(lun, err_alloc, 'meshProp', &
                 '  type mesh member surface (array)',res)    
         end if
         if (.not. allocated(this%normal)) then
            allocate(this%normal(3,this%nface),stat=res)
            if(res .ne. 0) rc = IOerr(lun, err_alloc, 'meshProp', &
                 '  type mesh member normal (array)',res)
         end if
         if (.not. allocated(this%bar_edge)) then
            allocate(this%bar_edge(3,this%nedge),stat=res)
            if(res .ne. 0) rc = IOerr(lun, err_alloc, 'meshProp', &
                 '  type mesh member bar_edge (array)',res)
         end if
      end if
      if (this%cnc_built .ge. 2) then
         ! nothing to be done here at this moment
      end if

      this%size_node = zero
      do icell=1,this%ncell
         n1=this%topol(1,icell)
         n2=this%topol(2,icell)
         n3=this%topol(3,icell)
         n4=this%topol(4,icell)
         s1n1=this%coord(1,n1)
         s2n1=this%coord(2,n1)
         s3n1=this%coord(3,n1)
         s1n2=this%coord(1,n2)
         s2n2=this%coord(2,n2)
         s3n2=this%coord(3,n2)
         s1n3=this%coord(1,n3)
         s2n3=this%coord(2,n3)
         s3n3=this%coord(3,n3)
         s1n4=this%coord(1,n4)
         s2n4=this%coord(2,n4)
         s3n4=this%coord(3,n4)
         this%bar_cell(1,icell)=(s1n1+s1n2+s1n3+s1n4)*onefourth
         this%bar_cell(2,icell)=(s2n1+s2n2+s2n3+s2n4)*onefourth
         this%bar_cell(3,icell)=(s3n1+s3n2+s3n3+s3n4)*onefourth
         this%size_cell(icell)=onesixth*abs(&
              &  (s1n4-s1n1)*((s2n2-s2n1)*(s3n3-s3n1)-(s3n2-s3n1)*(s2n3-s2n1))&
              & +(s2n4-s2n1)*((s3n2-s3n1)*(s1n3-s1n1)-(s1n2-s1n1)*(s3n3-s3n1))&
              & +(s3n4-s3n1)*((s1n2-s1n1)*(s2n3-s2n1)-(s2n2-s2n1)*(s1n3-s1n1)))

         ! area node
         do iloc=1,this%nnodeincell
            inode=this%topol(iloc,icell)
            this%size_node(inode) = this%size_node(inode) + &
                 this%size_cell(icell) * onefourth
         end do

      end do

      if (this%cnc_built .ge. 1) then

         do iface=1,this%nface_bc
            n1=this%face_iside(1,iface)
            n2=this%face_iside(2,iface)
            n3=this%face_iside(3,iface)
            s1n1=this%coord(1,n1)
            s2n1=this%coord(2,n1)
            s3n1=this%coord(3,n1)
            s1n2=this%coord(1,n2)
            s2n2=this%coord(2,n2)
            s3n2=this%coord(3,n2)
            s1n3=this%coord(1,n3)
            s2n3=this%coord(2,n3)
            s3n3=this%coord(3,n3)
            ! compute the distances between n1 and n2 (r)  and between n1 and n3 (d)
            r1=s1n1-s1n2
            r2=s2n1-s2n2
            r3=s3n1-s3n2
            d1=s1n1-s1n3
            d2=s2n1-s2n3
            d3=s3n1-s3n3

            this%surface_face(iface)=sqrt((r2*d3-r3*d2)*(r2*d3-r3*d2) &
                                         +(r3*d1-r1*d3)*(r3*d1-r1*d3) &
                                         +(r1*d2-r2*d1)*(r1*d2-r2*d1))*onehalf
            this%normal(1,iface)=(r2*d3-r3*d2)/this%surface_face(iface)
            this%normal(2,iface)=(-r1*d3+r3*d1)/this%surface_face(iface)
            this%normal(3,iface)=(r1*d2-r2*d1)/this%surface_face(iface)

!!$            if (proj_bar(this%normal(:,iedg), &
!!$               this%normal(:,iedg)=-this%normal(:,iedg)
!!$            end if
         end do

         do iface=this%nface_bc+1,this%nface
            n1=this%face_iside(1,iface)
            n2=this%face_iside(2,iface)
            n3=this%face_iside(3,iface)
            s1n1=this%coord(1,n1)
            s2n1=this%coord(2,n1)
            s3n1=this%coord(3,n1)
            s1n2=this%coord(1,n2)
            s2n2=this%coord(2,n2)
            s3n2=this%coord(3,n2)
            s1n3=this%coord(1,n3)
            s2n3=this%coord(2,n3)
            s3n3=this%coord(3,n3)
            ! compute the distances between n1 and n2 (r)  and between n1 and n3 (d)
            r1=s1n1-s1n2
            r2=s2n1-s2n2
            r3=s3n1-s3n2
            d1=s1n1-s1n3
            d2=s2n1-s2n3
            d3=s3n1-s3n3

            this%surface_face(iface)=sqrt((r2*d3-r3*d2)*(r2*d3-r3*d2) &
                                         +(r3*d1-r1*d3)*(r3*d1-r1*d3) &
                                         +(r1*d2-r2*d1)*(r1*d2-r2*d1))*onehalf
            this%normal(1,iface)=(r2*d3-r3*d2)/this%surface_face(iface)
            this%normal(2,iface)=(-r1*d3+r3*d1)/this%surface_face(iface)
            this%normal(3,iface)=(r1*d2-r2*d1)/this%surface_face(iface)

         end do

         do icell=1,this%ncell
            this%surface_tetra(icell)=zero
            do iloc=1,4
               iface=this%face_cnc(iloc,icell)
               this%surface_tetra(icell)=&
                    this%surface_tetra(icell)+this%surface_face(iface)
            end do
         end do

         do iedg=1,this%nedge
            n1=this%iside_edge(1,iedg)
            n2=this%iside_edge(2,iedg)
            s1n1=this%coord(1,n1)
            s2n1=this%coord(2,n1)
            s3n1=this%coord(3,n1)
            s1n2=this%coord(1,n2)
            s2n2=this%coord(2,n2)
            s3n2=this%coord(3,n2)
            this%bar_edge(1,iedg)=(s1n1+s1n2)*onehalf
            this%bar_edge(2,iedg)=(s2n1+s2n2)*onehalf
            this%bar_edge(3,iedg)=(s3n1+s3n2)*onehalf
         end do
            
      end if
    contains
      function proj_bar(normal,bar_edge,bar_cell) result(proj)
        implicit none
        ! lapack
        real(kind=double) :: ddot
        real(kind=double), intent(in  ) :: normal(2),bar_edge(2),bar_cell(2)
        real(kind=double)               :: proj
        proj=ddot(2,normal,1,bar_cell,1)-ddot(2,normal,1,bar_edge,1)
        
      end function proj_bar
    end subroutine meshProp

    !>---------------------------------------------------------------------
    !> procedure for averaging in a coarser mesh a quantity defined on a 
    !> refined mesh, weighting with respect to the volume of each
    !> refined tetrahedron, thus it works also if the mesh is not uniformly
    !> refined
    !>
    !>
    !> usage:    call  var%avg(grid, data, data_avg)
    !>
    !> where:
    !> \param[inout] subgrid  -> class(mesh3d), refined mesh
    !> \param[inout] grid     -> class(mesh3d), coarser mesh
    !> \param[in   ] data     -> real(ncell_sub) data to be averaged
    !> \param[out  ] data_avg -> real(ncell_sub) data averaged
    !>---------------------------------------------------------------------
    subroutine avg(subgrid, grid, data, data_avg)
      use Globals
      implicit none
      class(mesh3d),       intent(inout) :: subgrid, grid
      real(kind=double), intent(in   ) :: data(subgrid%ncell)
      real(kind=double), intent(inout) :: data_avg(grid%ncell)
      !local 
      integer :: icell,icell_sub
      data_avg=zero
      do icell_sub = 1, subgrid%ncell
         icell = subgrid%cell_parent(icell_sub)
         data_avg(icell) = data_avg(icell) + data(icell_sub)* &
              subgrid%size_cell(icell_sub)/grid%size_cell(icell)
      end do
    end subroutine avg

    !>------------------------------------------------------------------
    !> procedure for averaging in  coarser mesh a quantity defined on a 
    !> refined mesh.
    !>      $(data_avg(i)=sum_{i=1,4}data(subtriangles(i))/4$
    !> It works  only if triang_sub is the topology 
    !> of a conformally refined mesh
    !>
    !> usage:    call  avg_vec(subgrid, data, data_avg)
    !>
    !> where:
    !> \param[in ] ndata      -> integer. Number of data to be averaged
    !> \param[in ] data       -> real Dimension (ndata,ncell_sub) 
    !>                             data to be averaged
    !> \param[out] data_avg   -> real Dimension(ndata,ncell_sub) 
    !>                             data averaged
    !>------------------------------------------------------------------
    subroutine avg_vec(subgrid, grid,ndata,data, data_avg)
      use Globals
      implicit none
      class(mesh3d),       intent(in   ) :: subgrid
      class(mesh3d),       intent(in   ) :: grid
      integer,           intent(in   ) :: ndata
      real(kind=double), intent(in   ) :: data(ndata,subgrid%ncell)
      real(kind=double), intent(inout) :: data_avg(ndata,subgrid%ncell_parent)
      !local 
      integer :: icell,icell_sub

      data_avg=zero
      do icell_sub = 1, subgrid%ncell
         icell = subgrid%cell_parent(icell_sub)
         data_avg(:,icell) = data_avg(:,icell) + &
              data(:,icell_sub) * &
              subgrid%size_cell(icell_sub)/grid%size_cell(icell)
      end do


    end subroutine avg_vec


    !>------------------------------------------------------------------
    !> procedure for projection over a rifened mesh vreal var. data
    !> (public procedure for type mesh)
    !>
    !> usage:    call  avg(subgrid,data,data_avg)
    !>
    !> where:
    !> \param[in ] data       -> real(ncell_sub) data to be averaged
    !> \param[out] data_avg   -> real(ncell_sub) data averaged
    !>------------------------------------------------------------------
    subroutine proj_subgrid(subgrid, data, data_prj)
      use Globals
      implicit none
      class(mesh3d),       intent(in   ) :: subgrid
      real(kind=double), intent(in   ) :: data(subgrid%ncell/4)
      real(kind=double), intent(inout) :: data_prj(subgrid%ncell)
      !local 
      integer :: icell_sub

      do icell_sub = 1, subgrid%ncell 
         data_prj(icell_sub) = data(subgrid%cell_parent(icell_sub))
      end do

    end subroutine proj_subgrid

    !>------------------------------------------------------------------
    !> procedure for projection over a rifened mesh vreal var. data
    !> (public procedure for type mesh)
    !>
    !> usage:    call  avg(subgrid,data,data_avg)
    !>
    !> where:
    !> \param[in ] data       -> real(ncell_sub) data to be averaged
    !> \param[out] data_avg   -> real(ncell_sub) data averaged
    !>------------------------------------------------------------------
    subroutine projnode_subgrid(subgrid, data, data_prj)
      use Globals
      implicit none
      class(mesh3d),       intent(in   ) :: subgrid
      real(kind=double), intent(in   ) :: data(subgrid%nnode_parent)
      real(kind=double), intent(inout) :: data_prj(subgrid%nnode)
      !local 
      integer :: inode_sub

      do inode_sub = 1, subgrid%nnode
         data_prj(inode_sub) = onehalf * (&
              data(subgrid%node_parent(1,inode_sub)) + &
              data(subgrid%node_parent(2,inode_sub)))

      end do

    end subroutine projnode_subgrid

    !>------------------------------------------------------------------
    !> procedure 
    !> (public procedure for type mesh)
    !>
    !> usage:    call  
    !>
    !> where:
    !> \param[in ]   -> real(ncell_sub) data to be averaged
    !> \param[out]   -> real(ncell_sub) data averaged
    !>------------------------------------------------------------------
    subroutine heraldry(this, nend, ngrids, grids, tetra_ancestors)
      use Globals
      implicit none
      class(mesh3d),  intent(in   ) :: this
      integer   ,   intent(in   ) :: nend, ngrids
      type(mesh3d),   intent(in   ) :: grids(ngrids)
      integer,      intent(out) :: tetra_ancestors(grids(nend)%ncell)
      !local 
      integer :: icell, icell_father,grid_level, ncell_max

      ! max level of refined
      ncell_max = grids(nend)%ncell    
      do icell = 1, ncell_max
         tetra_ancestors(icell) = icell
      end do

      grid_level = grids(nend)%grid_level
      do while ( grid_level .ne. this%grid_level )
         do icell = 1, ncell_max
            icell_father = & 
                 grids(grid_level+1)%cell_parent(tetra_ancestors(icell))
            tetra_ancestors(icell) = icell_father 
         end do
         grid_level = grid_level - 1
      end do
    end subroutine heraldry

    !>------------------------------------------------------------------
    !> Evaluates the p-norm of functions piecewise constant on tetrahedrons
    !>------------------------------------------------------------------
    function normp_cell(this,power,data) result(total)
      use Globals
      implicit none
      class(mesh3d),       intent(in) :: this
      real(kind=double), intent(in) :: power
      real(kind=double), intent(in) :: data(this%ncell)
      real(kind=double) :: total
      !local 
      integer :: icell
      real(kind=double) :: ddot

      if ( power .eq.  0 ) then
         total = maxval( abs( data ) )
      else if ( power .eq. 1 ) then
         total = ddot(this%ncell, abs(data) ,1,this%size_cell,1)
      else
         total = zero
         do icell = 1, this%ncell
            total = total + abs(data(icell)) ** power * this%size_cell(icell)
         end do
         total = total ** ( one / power )
      end if
    end function normp_cell

    !>-----------------------------------------------------
    !> Find the index of the closer node to given point
    !>-----------------------------------------------------
    function closer_node(this,point) result(id_closernode)
      use Globals
      implicit none
      class(mesh3d), intent(in) :: this
      real(kind=double) ,intent(in) :: point(2)
      integer :: id_closernode
      !local 
      integer i
      real(kind=double) :: dist,dist_min
      real(kind=double) :: temp(2)
      dist_min=1.0d30
      do i = 1, this%nnode
         temp = point - this%coord(:,i)
         dist = sqrt(temp(1)**2+temp(2)**2)
         if ( dist < dist_min ) then
            id_closernode=i
            dist_min = dist
         end if
      end do
    end function closer_node


    !>-----------------------------------------------------
    !> Minimal mesh parameter
    !>-----------------------------------------------------
    function meshpar(this,flag) result(h)
      use Globals
      implicit none
      class(mesh3d), intent(in) :: this
      integer,     intent(in) :: flag 
      real(kind=double) :: h
      !local 
      if (flag .eq. 0) h=(6*sum(this%size_cell)/this%ncell)**onethird
      if (flag .eq. 1) h=minval((this%size_cell)**onethird/6.0d0)
      if (flag .eq. 2) h=maxval((this%size_cell)**onethird/6.0d0)
    end function meshpar
    !>-----------------------------------------------------
    !> Given a point gives the index of the closer
    !> point, the minimal_distance, the index of the supporting tetrahedron
    !> (index equal = 0 if dist_min < small ) 
    !>-----------------------------------------------------
    subroutine info_point(this,point,&
         index_closer_node, index_closer_tria, dist_min)
      use Globals
      implicit none
      class(mesh3d),        intent(in ) :: this
      real(kind=double) , intent(in ) :: point(2)
      integer,            intent(out) :: index_closer_node
      integer,            intent(out) :: index_closer_tria
      real(kind=double) , intent(out) :: dist_min
      !local 
      integer j,inode,icell
      real(kind=double) :: dist
      real(kind=double) :: temp(2)
      dist_min=1.0d30
      do icell = 1, this%ncell
         do j=1,3
            inode = this%topol(j,icell)
            temp = point - this%coord(:,inode)
            dist = sqrt(temp(1)**2+temp(2)**2)
            if ( dist < dist_min ) then
               index_closer_node=inode
               index_closer_tria=icell
               dist_min = dist
               if (dist_min < small) index_closer_tria=0
            end if
         end do
      end do
    end subroutine info_point
    !>-------------------------------------------------------------------
    !> Procedure that reorders mesh nodes
    !> given a permutation vector.
    !> It works on members triang and coord
    !> and rebuilds all other member if necessary. 
    !>
    !> usage:
    !>     call 'var'%reorder(lun_err,perm_nodes,inv_perm_nodes)
    !>
    !> where:
    !> \param[in] lun_err        -> integer. I/O unit for error message output
    !> \param[in] perm_nodes     -> integer. dimension(nnode)
    !>                                       node permutation array
    !> \param[in] inv_perm_nodes -> integer. dimension(nnode)
    !>                                       inverse of perm_nodes
    !>------------------------------------------------------------------
    subroutine reorder(this,lun_err,perm_nodes,inv_perm_nodes)
      use Globals
      implicit none
      class(mesh3d), intent(inout) :: this
      integer,     intent(in   ) :: lun_err
      integer,     intent(in   ) :: perm_nodes(this%nnode) 
      integer,     intent(in   ) :: inv_perm_nodes(this%nnode) 
      ! local 
      logical :: rc
      integer :: res
      integer :: inode,icell,i,iedge
      integer :: ncell, nnodeincell
      type(mesh3d) :: local_mesh

      ncell       = this%ncell
      nnodeincell = this%nnodeincell

      !  Permute the nodes according to the permutation vector.
      call double_col_permute(3, this%nnode, perm_nodes, this%coord)
            
      !
      !  Permute the node indices in the tetrahedron array.
      do icell = 1, ncell
         do i = 1, nnodeincell
            inode = this%topol(i,icell)
            this%topol(i,icell) = inv_perm_nodes ( inode )
         end do
      end do

      ! By UNANIMOUS DECISION we check in connection_faces
      ! if arrays are allocated

      ! Create local working mesh
      select case (this%cnc_built)
      case (0)
         ! Rebuild real Properties
         call this%meshProp(lun_err)
      case (1)
         ! Rebuild 1st level connections
         call this%build_face_connection(lun_err)
         ! Rebuild real Properties
         call this%meshProp(lun_err)
      end select

      ! Renumber also node_parent if this is a subgrid
      if ( this%grid_level .gt. 0 ) then
         call integer_col_permute(2,this%nnode,perm_nodes, &
              this%node_parent)
      end if

    end subroutine reorder
   
    !>-------------------------------------------------------------------
  !> Procedure to build the permutation of the nodes to
  !> minimize the matrix bandwidth for Galerkin P1 FEM
  !>
  !> usage:
  !>     call 'var'%build_rcm_permutation(lun_err,
  !>                                      perm_nodes,
  !>                                      inv_perm_nodes)
  !>
  !> where:
  !> \param[in] lun_err        -> integer. I/O unit for error message output
  !> \param[in] perm_nodes     -> integer. dimension(nnode)
  !>                                       node permutation array
  !> \param[in] inv_perm_nodes -> integer. dimension(nnode)
  !>                                       inverse of perm_nodes
  !>-------------------------------------------------------------------
  subroutine build_rcm_permutation(this,lun_err,perm_nodes,iperm_nodes)
    use Globals
    use SparseMatrix
    implicit none
    class(mesh3d),        intent(inout) :: this
    integer,            intent(in   ) :: lun_err
    integer, optional,  intent(inout) :: perm_nodes(this%nnode) 
    integer, optional,  intent(inout) :: iperm_nodes(this%nnode) 
    !local 
    logical :: rc
    integer :: res,i,j,inode
    type(spmat) :: adj_matrix
    
    perm_nodes=0
    iperm_nodes=0

    call this%nodenode_connection(lun_err,.False.,adj_matrix) 
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  ADJ bandwidth = ', adj_matrix%bandwidth()
    !
    !  Compute the RCM permutation.
    !
    call adj_matrix%genrcm ( lun_err, perm_nodes,iperm_nodes )
    !
    !write(*,*) iperm_nodes
    
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Permuted ADJ bandwidth = ', &
         adj_matrix%bandwidth(perm_nodes,iperm_nodes)
    write ( *, '(a)' ) ' '
    
   
    !
    !  Free memory.
    !
    call adj_matrix%kill(lun_err)   
  end subroutine build_rcm_permutation

  subroutine nodenode_connection(this, lun_err,include_diagonal,adj_matrix)
    use Globals
    use SparseMatrix
    implicit none
    class(mesh3d),  intent(inout) :: this
    integer,      intent(in   ) :: lun_err
    logical,      intent(in   ) :: include_diagonal
    type(spmat),  intent(inout) :: adj_matrix
    !local
    logical :: rc
    logical :: clean_memory
    integer :: res
    integer :: iedge,inode,i1,i2,nnode,nterm,j
    integer, allocatable :: iaw(:)
    integer, allocatable :: jaw(:)
    integer, allocatable :: count(:)
    
    if ( .not. allocated(this%iside_edge))  then
       call this%connection_edges(lun_err)
       clean_memory = .true.
    end if

    nnode = this%nnode
    if ( include_diagonal) then
       nterm = 2*this%nedge+this%nnode
    else
       nterm = 2*this%nedge
    end if

    allocate(iaw(nnode+1),jaw(nterm),count(nnode),stat=res)
    if(res .ne. 0) rc = IOerr(lun_err, err_alloc, 'init_star', &
         'work array ia ja count',res)
    
    ! build ia pointer
    if ( include_diagonal) then
       count=1
    else
       count=0
    end if
    do iedge = 1,this%nedge
       i1=this%iside_edge(1,iedge)
       i2=this%iside_edge(2,iedge)
       count(i1) = count(i1)+1
       count(i2) = count(i2)+1
       iaw(this%iside_edge(1,iedge)) = iaw(this%iside_edge(1,iedge))+1 
       iaw(this%iside_edge(2,iedge)) = iaw(this%iside_edge(2,iedge))+1
    end do
    iaw(1) = 1 
    do inode = 1, nnode
       iaw(inode+1) =  iaw(inode) + count(inode) 
    end do
       
    if ( include_diagonal ) then
       ! add self connection
       do inode = 1, this%nnode
          jaw(iaw(inode)) = inode
       end do
       count = 1
    else
       count = 0
    end if
    
    do iedge = 1, this%nedge
       i1 = this%iside_edge(1,iedge)
       i2 = this%iside_edge(2,iedge)
       jaw( iaw(i1)+count(i1) ) = i2 
       count( i1 ) = count( i1 ) + 1 
       jaw( iaw(i2)+count(i2) ) = i1
       count( i2 ) = count( i2 ) + 1 
    end do

   
    call adj_matrix%init(lun_err,&
         nnode,nnode, nterm,&
         storage_system='csr',&
         is_symmetric=.true.)

    adj_matrix%is_sorted = .False.    
    adj_matrix%ia    = iaw
    adj_matrix%ja    = jaw
    adj_matrix%coeff = zero
   
    call adj_matrix%sort()


    ! add extra diagonal terms
    deallocate(iaw,jaw,count,stat=res)
    if(res .ne. 0) rc = IOerr(lun_err, err_dealloc, 'init_star', &
         'work array ia ja count',res)
    
    if (clean_memory) then
       ! clean array built in connection
       deallocate(this%edge_cnc,this%iside_edge,stat=res)
       if (res.ne.0) rc=IOerr(lun_err, err_dealloc, 'kill_mesh', &
            'mesh members edge_cnc this%iside_edge')
    end if

  end subroutine nodenode_connection

   !>-------------------------------------------------
  !> Precedure for computing 
  !> lapcian(r) = \int_{T_r} (\Delta u ) v_r
  !> with 
  !>  u   : piecewise constant fucntions of the triangulation
  !>  v_r : characteristion function of element r
  !>
  !> ( public procedure for type tdpotsys)
  !> 
  !> usage: call var%build_nrm_grad_vars(ctrl)
  !>
  !<-------------------------------------------------
  subroutine p0_laplacian(this,pot,laplacian)
    use Globals
    implicit none
    class(mesh3d), intent(in) :: this
    real(kind=double),intent(in) :: pot(this%ncell)
    real(kind=double),intent(out) :: laplacian(this%ncell)
    !local
    integer :: iface,icell1,icell2
    real(kind=double) :: vec(3),diff,dist_bar
    real(kind=double) :: dnrm2
    
    !
    ! internal edges
    !
    laplacian = zero
    do iface = 1, this%nface - this%nface_bc
       icell1 = this%face_iside(1,iface)
       icell2 = this%face_iside(2,iface)
       vec = this%bar_cell(:,icell1)-this%bar_cell(:,icell2)
       dist_bar=dnrm2(3,vec,1)
       diff = pot(icell1) - pot(icell2)
       laplacian(icell1) = laplacian(icell1) + &
            diff / dist_bar * this%surface_face(iface)  
       laplacian(icell2) = laplacian(icell2) + &
            diff / dist_bar * this%surface_face(iface)  
    end do
    
    !
    ! no boudary contribution
    !
    
    

  end subroutine p0_laplacian


end module Geometry3d
