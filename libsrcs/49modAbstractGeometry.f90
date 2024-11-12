module AbstractGeometry
  use Globals
  implicit none
  private
  public :: abs_simplex_mesh, mesh_constructor_from_data, mesh_destructor, mesh_constructor_from_cell_selection, pure_cartesian2barycentric
  type :: abs_simplex_mesh
     !> Logical flag if type has been initialized
     logical :: is_initialized=.False.
     !> Dimension of the space where grid is embedded
     !> ambient_dimension=2 flat 2d grid (coord(3,:) = 0)
     !> ambient_dimension=3 3d grid or surface grid
     integer :: ambient_dimension=0
     !> Mesh dimension as manifold  
     !> logical_dimension=1 cell=edge
     !> logical_dimension=2 cell=triangle
     !> logical_dimension=3 cell=tetrahedron
     integer :: logical_dimension=0
     !> Label for grid
     !> 'edge'        = 1d cell
     !> 'triangle'    = triangles ( with 2d or 3d coordinates)
     !> 'tetrahedron' = tetrahedron 
     !> 'square'      = flat square aligned with xy axis
     character(len=256) :: cell_type='empty'
     !> Cell type as vtk format
     !> 3  = edge
     !> 5  = triangles 
     !> 10 = tetrahedron
     !> 9  = square
     integer :: cell_id=0    
     !> Label for grid
     !> '2d'      = 2d flat grid  composed by triangles  coord(3,:)=0.0
     !> 'surface' = 2d grid embedded in 3d composed by triangles
     !> '3d'      = 3d grid composed by tetrahedrons
     character(len=256) :: mesh_type
     !> Label for the normal stored
     !> normal_storage='node' : normal_node is allocated and used
     !> normal_storage='cell' : normal_cell is allocated and used
     character(len=4) :: normal_storage='cell'
     !> Number of nodes
     integer :: nnode=0
     !> Number of cells
     integer :: ncell=0
     !> Number of node for each 
     integer :: nnodeincell=0
     !> Number of zones
     integer :: nzone=0
     !> Dimension (3,nnode). 
     !> Nodal coordinate. 
     real(kind=double), allocatable :: coord(:,:)
     !> Topological info
     !> Dimension (nnode_in_cell,ncell).
     integer, allocatable :: topol(:,:) 
     !> Dimension (ncell).
     !> Size of the cell
     !> **Computed in procedure Prop.**
     real(kind=double), allocatable :: size_cell(:)
     !> Dimension (nnode).
     !> Nodal area/volume
     !> **Computed in procedure Prop.**
     real(kind=double), allocatable :: size_node(:)
     !> Dimension (3,ncell).
     !> Normal to the cell
     !> **Computed in procedure Prop.**
     real(kind=double), allocatable :: normal_cell(:,:)
     !> Dimension (ncell).
     !> coordinates of center of mass of the cell
     real(kind=double), allocatable :: bar_cell(:,:)
     !>-------------------------------------------------
     !> Edges quantities
     !>-------------------------------------------------
     !> Number of edges
     integer :: nedge=0
     !> Number of boundary edges
     integer :: nedge_bc=0
     !> Dimension (nnodedgeincell,ncell). 
     !> Edges (in the global enumeration) in each cell.
     !> **Calculated in the constructor.**
     integer, allocatable :: edge_cnc(:,:)
     !> Dimension (2,nedge). 
     !> Nodes for each edge. 
     !> **Calculated in the constructor.**
     integer, allocatable :: edge_iside(:,:)
     !> Dimension (2,nedge). 
     !> Left-Right cell for each edge. 
     integer, allocatable :: edge_plist(:,:)
     !> Dimension (3,ncell). 
     !> Neighborhood cells for each cells 
     integer, allocatable :: neigh_cell(:,:)
     !> Dimension (nedge). 
     !> Length of 1d-edges 
     integer, allocatable :: size_edge(:)
     !> Number of boundary nodes
     integer :: nnode_bc=0
     !> Dimension (nnode_bc). 
     !> List of boundary nodes 
     integer, allocatable :: node_bc(:)
     !> Number of boundary cell
     integer :: nboundary_cell=0
     !> Dimension (nboundary_cell). 
     !> List of boundary nodes 
     integer, allocatable :: boundary_cell(:)
     !>-----------------------------------------------
     !> Faces quantities 
     !>-----------------------------------------------
     !> Number of mesh faces (including boundary). 
     !> Calculated in the constructor.
     integer :: nface = 0
     !> Number of boundary faces.
     !> **Calculated in the constructor.**
     integer :: nface_bc = 0
     !> Dimension (4,ncell). 
     !> Faces (in the global enumeration) in each cell.
     !> **Calculated in the constructor.**
     integer, allocatable :: face_cnc(:,:)
     !> Dimension (4,ncell).
     !> Neighboring tetrahedrons for each tetrahedron.
     !> **Calculated in the constructor.**
     integer, allocatable :: neigh_face(:,:)
     !> Dimension (3,nface).
     !> Nodes for each face. 
     !> **Calculated in the constructor.**
     integer, allocatable :: face_iside(:,:)
     !> Dimension (2,nface). 
     !> East and West cell for each face
     !> East and West are calculated with respect to the nodal ordering in CELL.
     !> **Calculated in the constructor.**
     integer, allocatable :: face_plist(:,:)
     !>-----------------------------------------------
     !> Subgrid quantities 
     !>-----------------------------------------------
     !> Grid level 
     integer  :: grid_level=0
     !> Number of nodes of parent grid
     integer :: nnode_parent=0
     !> Number of tetrahedrons parent grid = ncell/8
     integer :: ncell_parent=0
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
     integer, allocatable :: node_parent(:,:)
     !> Dimension (ncell)
     !> Parent tetrahedrons
     integer, allocatable :: cell_parent(:)
   contains
     !> Write basic property of mesh
     !> (procedure public for type abs_simplex_mesh)
     procedure, public, pass :: info => info_mesh
     !> Free memory procedure
     !> (procedure public for type abs_simplex_mesh)
     procedure, public, pass :: kill => kill_mesh 
     !> Read mesh from file
     !> (procedure public for type abs_simplex_mesh)
     procedure, public, pass :: read_mesh
     !> Write mesh to file
     !> (procedure public for type abs_simplex_mesh)
     procedure, public, pass :: write_mesh
     !> Build a smaller grid removing list of nodes or cells
     !> (procedure public for type abs_simplex_mesh)
     procedure, public, pass :: init_selection
     !> Build a smaller grid removing list of nodes or cells
     !> (procedure public for type abs_simplex_mesh)
     procedure, public, pass :: init_from_data
     !> Define the sparse (nell x nnode) matrix 
     !> wiht the node connecctionity
     !> (procedure public for type abs_simplex_mesh)
     procedure, public, pass :: build_connection_matrix
     !> Define the sparse (nnode x nnode) matrix 
     !> wiht the node-node connectiovity
     !> (procedure public for type abs_simplex_mesh)
     procedure, public, pass :: build_nodenode_matrix
     !> Define the sparse (ncell x ncell) matrix
     !> with the cell-cell connectiovity
     !> (procedure public for type abs_simplex_mesh)
     procedure, public, pass :: build_cellcell_matrix
     !> Define the sparse (ncell x ncell) matrix
     !> with the cell-cell connectiovity
     !> (procedure public for type abs_simplex_mesh)
     procedure, public, pass :: build_cellcell_graph
     !> Build size cell
     !> (procedure public for type abs_simplex_mesh)
     procedure, public, pass :: build_size_cell
     !> Build nodal area/volume
     !> (procedure public for type abs_simplex_mesh)
     procedure, public, pass :: build_size_node
     !> Build normal vector to cell
     !> (procedure public for type abs_simplex_mesh)
     procedure, public, pass :: build_normal_cell
     !> Build array cells barycenters
     !> (procedure public for type abs_simplex_mesh)
     procedure, public, pass :: build_bar_cell
     !> calculate Lebesgue p-norm of piecewise-constant array
     !> (public procedure for type mesh)
     procedure, public, pass :: normp_cell
     !> Build edges conenctions
     !> (procedure private for type abs_simplex_mesh)
     procedure, public, pass :: build_edge_connection
     !> Build faces conenctions
     !> (procedure private for type abs_simplex_mesh)
     procedure, public, pass :: build_face_connection
     !> Set number and list of boundary nodes.
     !> (procedure public for type abs_simplex_mesh)
     procedure, public, pass :: build_nodebc
     !> Set number and list of boundary nodes.
     !> (procedure public for type abs_simplex_mesh)
     procedure, public, pass :: build_boundary_cell
     !> Build edges conenctions
     !> (procedure private for type abs_simplex_mesh)
     procedure, public, pass :: renumber => renumber_mesh
     !>-----------------------------------------------
     !> Subgrid procedures
     !>-----------------------------------------------
     !> Interpolate data defined on cells on a coaser 
     !> grid to data defined on cells on finer grids
     !> (public procedure for type mesh)
     procedure, public, pass :: proj_subgrid
     !> Interpolate data defined on nodes on a coaser 
     !> grid to data defined on nodes on finer grids
     !> (public procedure for type mesh)
     procedure, public, pass :: projnode_subgrid
     !> Project data defined on cell of a coaser grid 
     !> to data defined on cells of a finer grid
     !> (public procedure for type mesh)
     procedure, public, pass :: avg_cell_subgrid
     !> Procedure to read subgrid cell_parent and node_parent
     !> (public procedure for type mesh)
     procedure, public, pass :: read_parent => read_grid_subgrid_relation
     !> Procedure to write subgrid cell_parent and node_parent
     !> (public procedure for type mesh)
     procedure, public, pass :: write_parent => write_grid_subgrid_relation
     !> builds mesh refining a previous mesh 
     !> (private procedure for type mesh)
     procedure, public, pass :: refine => build_refined_topol_coord
     !> Procedure to reorder nodes to reduce bandwidth
     !> Reverse Cuthill-Mckee 
     !> (private procedure for type mesh)
     procedure, public, pass :: RCM_node_reorder
     !>-----------------------------------------------
     !> node2cell and cell2node
     !>-----------------------------------------------
     !> Procedure for computing the cell values of a function
     !> given the nodal values
     !> (procedure public for type abs_simplex_mesh)
     procedure, public, pass :: eval_node2cell
     !> Procedure for computing the nodal values of a function
     !> given the cell values
     !> (procedure public for type abs_simplex_mesh)
     procedure, public, pass :: eval_cell2node
     !> Procedure for computing the nodal values of a function
     !> given the material zone values
     !> (procedure public for type abs_simplex_mesh)
     procedure, public, pass :: eval_zone2node
     !> Procedure to compute different tipical mesh parameter
     !> (procedure private for type abs_simplex_mesh)
     procedure, public, pass :: meshpar     
     !> Write basic property of mesh
     !> (procedure public for type abs_simplex_mesh)
     procedure, public, pass :: scalar_product
     !> Return the list of the nodes in the cell given a cell index
     !> (procedure private for type abs_simplex_mesh)
     procedure, public, pass :: get_coord_of_node
     !> Return the list of the nodes in the cell given a cell index
     !> (procedure private for type abs_simplex_mesh)
     procedure, public, pass :: get_nnodes_in_cell
     !> Return the list of the nodes in the cell given a cell index
     !> (procedure private for type abs_simplex_mesh)
     procedure, public, pass :: get_nodes_in_cell
     !> Build transformation from cartesian to barycentry coordinate
     !> (procedure public for type abs_simplex_mesh)
     procedure, public, pass :: cartesian2barycentric
  end type abs_simplex_mesh


  type, extends(abs_simplex_mesh) :: metric_mesh
     !> Three vector
     real(kind=double) :: metric(2,2)
   contains
     !> Write basic property of mesh
     !> (procedure public for type abs_simplex_mesh)
     procedure, public, pass :: scalar_product => metric_scalar_product

  end type metric_mesh
contains
  !>-------------------------------------------------------------
  !> Print basic grid informations.
  !> (procedure private for type abs_simplex_mesh)
  !>
  !> usage:
  !>     call 'varr'%info(lun)
  !>
  !> where:
  !> \param[in] lun  -> integer. IO unit 
  !<-------------------------------------------------------------
  subroutine info_mesh(this, lun)
    use Globals
    implicit none
    !vars
    class(abs_simplex_mesh),    intent(in) :: this
    integer,        intent(in ) :: lun

    write(lun,*) ' Number of nodes         = ', this%nnode
    write(lun,*) ' Number of cells         = ', this%ncell
    write(lun,*) ' Cell type               = ', etb(this%cell_type)
    write(lun,*) ' Ambient dimension       = ', this%ambient_dimension
    write(lun,*) ' Logical dimension       = ', this%logical_dimension

  end subroutine info_mesh

  !>-------------------------------------------------------------
  !> Part of the static constructor.
  !> (procedure private for type abs_simplex_mesh)
  !> Instantiate and initilize reading from input file
  !> variable of type abs_simplex_mesh
  !>
  !> usage:
  !>     call this%read_mesh(lun_err, file2read)
  !>
  !> where:
  !> \param[in] lun_err   -> integer. IO unit for error msg
  !> \param[in] file2read -> type(file). I/O file information
  !<-------------------------------------------------------------
  subroutine read_mesh(this, lun_err, file2read)
    use Globals
    implicit none
    !vars
    class(abs_simplex_mesh),    intent(inout) :: this
    integer,        intent(in ) :: lun_err
    type(file),     intent(in ) :: file2read
    ! local vars
    integer :: u_number
    integer :: j, k, res, ncoord,nnodeincell,cell_id
    logical :: rc
    character(len=256) :: rdwr,str,fname,first_line,inputs
    integer :: ambient_dimension
    integer :: cells(5),cells_bis(5)
    real(kind=double) :: zmax

    !
    ! free memory 
    !
    if ( this%is_initialized ) call this%kill(lun_err)
    
    u_number    = file2read%lun
    fname       = file2read%fn

    read(u_number,*,iostat=res) this%nnode
    if(res .ne. 0) rc = IOerr(lun_err, err_inp , 'read_mesh', &
         etb(fname) // ' vars nnode',res)

    !
    ! Read second line. 
    ! First try reading two integers.
    ! If it fails read just the one integer and nnodeincell 
    ! is evaluted later.
    ! If it fails again file is corrupted
    !
    read(u_number,'(a)',iostat=res) inputs
    if(res .ne. 0) rc = IOerr(lun_err, err_inp , 'read_mesh', &
         etb(fname) // ' second line',res)

    this%cell_id=0
    read(inputs,*,iostat=res) this%ncell, nnodeincell, cell_id
    if ( res .eq. 0) then
       this%nnodeincell = nnodeincell
       this%cell_id     = cell_id
    else
       read(inputs,*,iostat=res) this%ncell, nnodeincell
       if ( res .eq. 0) then
          this%nnodeincell = nnodeincell 
       else
          nnodeincell=0
       end if
    end if


    allocate(this%coord(3,this%nnode),stat=res)
    if(res .ne. 0) rc = IOerr(lun_err, err_alloc, 'read_mesh', &
         '  type mesh member coord (array)',res)

    !
    ! read first line to check if the the third coordinate is present 
    ! 
    read(u_number,'(a)',iostat=res) first_line
    ambient_dimension=3
    zmax=0
    read(first_line,*,iostat=res) (this%coord(k,1),k=1,ambient_dimension)
    if (res .ne. 0) THEN
       this%coord(3,:) = zero
       ! read only first two column, the first is initialized to zero
       ambient_dimension = 2
       read(first_line,*,iostat=res) (this%coord(k,1),k=1,ambient_dimension)
       if(res .ne. 0) THEN
          write(rdwr,'(i5)') 1
          str=etb(rdwr)//'/'
          rc = IOerr(lun_err, err_inp , 'read_mesh', &
               trim(etb(fname)) //&
               ' type mesh member array coord at line '//&
               trim(str),res)
       end if
    end if
    zmax=max(zmax,this%coord(3,1))
    !
    ! read remaing lines
    !
    do j=2,this%nnode
       read(u_number,*,iostat=res) (this%coord(k,j),k=1,ambient_dimension)
       if(res .ne. 0) THEN
          write(rdwr,'(i5)') j
          str=etb(rdwr)//'/'
          rc = IOerr(lun_err, err_inp , 'read_mesh', &
               trim(etb(fname)) //&
               ' type mesh member array coord at line '//&
               trim(str),res)
       end if
       zmax=max(zmax,this%coord(3,j))
    end do


    !
    ! if the 3 componet is equal to zero change ambient dimension
    !
    if (abs(zmax) .lt. small ) ambient_dimension=2


    ! if nnodeincell is defined allocate and read 
    ! topol, otherwise deduced from number of integer in 
    ! first line
    ! read first line to check if the the third coordinate is present 
    ! 
    if ( nnodeincell .eq. 0 ) then
       read(u_number,'(a)',iostat=res) first_line
       if(res .ne. 0) rc = IOerr(lun_err, err_inp, 'read_mesh', &
            ' reading first line of topol to test nncell',res)
       res = -1
       !
       ! start from tetrahedron and go back
       !
       nnodeincell = 4
       do while ( (res .ne. 0 ) .and. (nnodeincell.ge.3))
          read(first_line,*,iostat=res) (cells(k), k=1,nnodeincell+1)
          if(res.ne. 0) nnodeincell=nnodeincell-1
       end do
       if ( nnodeincell .lt. 3 ) &
            rc = IOerr(lun_err, err_inp, 'read_mesh', &
            '  nnodeincell not set ',res)
       select case ( nnodeincell)
       case (4)
          this%cell_id = 10
       case (3)
          this%cell_id = 5
       case (2)
          this%cell_id = 3
       end select


       !
       ! copy first cell and than read the others
       !
       this%nnodeincell=nnodeincell
       allocate(this%topol(nnodeincell+1,this%ncell),stat=res)
       if(res .ne. 0) rc = IOerr(lun_err, err_alloc, 'read_mesh', &
            '  type mesh member triangles (array)',res)
       this%topol(1:nnodeincell+1,1)=cells(1:nnodeincell+1)
       do j=2,this%ncell
          read(u_number,*,iostat=res) (this%topol(k,j), k=1,nnodeincell+1)
          if(res .ne. 0) then
             write(rdwr,'(i7)') j
             str=etb(rdwr)//'/'
             write(rdwr,'(i7)') k
             str=trim(str)//etb(rdwr)
             rc = IOerr(lun_err, err_inp , 'read_mesh', &
                  trim(etb(fname)) // &
                  ' type mesh member array triang at line/col '//trim(str),res)
          end if
       end do
    else
       this%nnodeincell=nnodeincell
       allocate(this%topol(nnodeincell+1,this%ncell),stat=res)
       if(res .ne. 0) rc = IOerr(lun_err, err_alloc, 'read_mesh', &
            '  type mesh member triangles (array)',res)
       do j=1,this%ncell
          read(u_number,*,iostat=res) (this%topol(k,j), k=1,nnodeincell+1)
          if(res .ne. 0) then
             write(rdwr,'(i7)') j
             str=etb(rdwr)//'/'
             write(rdwr,'(i7)') k
             str=trim(str)//etb(rdwr)
             rc = IOerr(lun_err, err_inp , 'read_mesh', &
                  trim(etb(fname)) // &
                  ' type mesh member array triang at line/col '//trim(str),res)
          end if
       end do
    end if

    this%nzone = maxval(this%topol(nnodeincell+1,:),this%ncell)

    !
    ! set mesh characteristic
    !
    this%ambient_dimension=ambient_dimension
    select case (this%nnodeincell) 
    case (4) 
       select case (this%cell_id)
       case (9)
          this%cell_type='square'
       case default
          this%cell_type='tetrahedron'
       end select
       this%mesh_type='3d'
       this%logical_dimension = 3
    case (3) 
       this%cell_type='triangle'
       select case ( ambient_dimension )
       case (2) 
          this%mesh_type='2d'
          this%logical_dimension = 2
       case (3)
          this%mesh_type='surface'
          this%logical_dimension = 3
       end select
    case (2)
       this%cell_type='edge'
       this%mesh_type='graph'
       this%logical_dimension = 1
    end select

    this%is_initialized = .True.
  end subroutine read_mesh

  !>-------------------------------------------------------------
  !> Part of the static constructor.
  !> (procedure public for type abs_simplex_mesh)
  !> Write mesh into file
  !>
  !> usage:
  !>     call this%read_mesh(lun_err, file2write)
  !>
  !> where:
  !> \param[in] lun_err   -> integer. IO unit for error msg
  !> \param[in] file2write -> type(file). I/O file information
  !<-------------------------------------------------------------
  subroutine write_mesh(this, lun_err, file2write,out_format)
    use Globals
    implicit none
    !vars
    class(abs_simplex_mesh), intent(in) :: this
    integer,        intent(in ) :: lun_err
    type(file),     intent(in ) :: file2write
    character(len=*), optional, intent(in):: out_format 
    ! local vars
    integer :: u_number
    integer :: j, k, res, ncoord,nnodeincell
    logical :: rc
    character(len=256) :: rdwr,str,fname,first_line,inputs,str1,str2
    integer :: ambient_dimension
    integer :: cells(5)
    real(kind=double) :: zmax
    character(len=256) :: useformat

    u_number    = file2write%lun
    fname       = file2write%fn

    if ( present(out_format) )then
       useformat=(etb(out_format))
    else
       useformat='ascii'
    end if

    select case (useformat)
    case('ascii')
       !
       ! write head
       !
       write(u_number,*,iostat=res) this%nnode, '             !  # of nodes '
       if(res .ne. 0) rc = IOerr(lun_err, err_out , 'read_mesh', &
            etb(fname) // ' vars nnode',res)
       write(u_number,*,iostat=res) this%ncell, this%nnodeincell, this%cell_id, &
            ' ! # of cells, # of nodes in cell, #cell id ',& 
            ' mesh_type =', etb(this%mesh_type)
       if(res .ne. 0) rc = IOerr(lun_err, err_out , 'read_mesh', &
            etb(fname) // ' vars this%ncell, nnodeincell, mesh_type',res)

       !
       ! write coordinates
       !
       j=1
       write(u_number,*,iostat=res) (this%coord(k,j),k=1,3), ' ! x,y,z '
       if(res .ne. 0) THEN
          write(rdwr,'(i5)') j
          str=etb(rdwr)//'/'
          rc = IOerr(lun_err, err_out , 'read_mesh', &
               trim(etb(fname)) //&
               ' type mesh member array coord at line '//&
               trim(str),res)
       end if
       do j=2,this%nnode
          write(u_number,*,iostat=res) (this%coord(k,j),k=1,3)
          if(res .ne. 0) THEN
             write(rdwr,'(i5)') j
             str=etb(rdwr)//'/'
             rc = IOerr(lun_err, err_out , 'read_mesh', &
                  trim(etb(fname)) //&
                  ' type mesh member array coord at line '//&
                  trim(str),res)
          end if
       end do

       str1=''
       do j=1,this%nnodeincell
          write(str2,'(a,I1,a)') ' n',j,','
          str1=(etb(str1)//etb(str2))
       end do
       str1=(etb(str1)// ' material')
       !
       ! write topol
       !
       j=1
       write(u_number,*,iostat=res) (this%topol(k,j), k=1,this%nnodeincell+1), '    !',  etb(str1)
       if(res .ne. 0) then
          write(rdwr,'(i5)') j
          str=etb(rdwr)//'/'
          write(rdwr,'(i5)') k
          str=trim(str)//etb(rdwr)
          rc = IOerr(lun_err, err_inp , 'write_mesh', &
               trim(etb(fname)) // &
               ' type mesh member array triang at line/col '//trim(str),res)
       end if
       do j=2,this%ncell
          write(u_number,*,iostat=res) (this%topol(k,j), k=1,this%nnodeincell+1)
          if(res .ne. 0) then
             write(rdwr,'(i5)') j
             str=etb(rdwr)//'/'
             write(rdwr,'(i5)') k
             str=trim(str)//etb(rdwr)
             rc = IOerr(lun_err, err_inp , 'write_mesh', &
                  trim(etb(fname)) // &
                  ' type mesh member array triang at line/col '//trim(str),res)
          end if
       end do

    case ('gid')
       ! write head defining mesh type
       select case (this%nnodeincell)
       case(3)
          write(str,*) ' MESH dimension 3 ElemType Triangle Nnode 3'
       case (4)
          write(str,*) ' MESH dimension 3 ElemType Tetrahedron Nnode 4'
       end select
       write(u_number,*,iostat=res) etb(str)
       if(res .ne. 0) rc = IOerr(lun_err, err_out , 'read_mesh', &
            etb(fname) // ' vars nnode',res)

       write(u_number,*,iostat=res) 'Coordinates '
       if(res .ne. 0) rc = IOerr(lun_err, err_out , 'read_mesh', &
            etb(fname) // ' vars nnode',res)
       do j=1,this%nnode
          write(u_number,*,iostat=res) j, (this%coord(k,j),k=1,3)
          if(res .ne. 0) THEN
             write(rdwr,'(i5)') j
             str=etb(rdwr)//'/'
             rc = IOerr(lun_err, err_out , 'read_mesh', &
                  trim(etb(fname)) //&
                  ' type mesh member array coord at line '//&
                  trim(str),res)
          end if
       end do

       write(u_number,*,iostat=res) 'End Coordinates '
       if(res .ne. 0) rc = IOerr(lun_err, err_out , 'read_mesh', &
            etb(fname),res)

       write(u_number,*,iostat=res) ' '
       if(res .ne. 0) rc = IOerr(lun_err, err_out , 'read_mesh', &
            etb(fname),res)
       write(u_number,*,iostat=res) 'Elements '
       if(res .ne. 0) rc = IOerr(lun_err, err_out , 'read_mesh', &
            etb(fname),res)
       do j=1,this%ncell
          write(u_number,*,iostat=res) j, (this%topol(k,j), k=1,this%nnodeincell)
          if(res .ne. 0) then
             write(rdwr,'(i5)') j
             str=etb(rdwr)//'/'
             write(rdwr,'(i5)') k
             str=trim(str)//etb(rdwr)
             rc = IOerr(lun_err, err_inp , 'write_mesh', &
                  trim(etb(fname)) // &
                  ' type mesh member array triang at line/col '//trim(str),res)
          end if
       end do


       write(u_number,*,iostat=res) 'End Elements '
       if(res .ne. 0) rc = IOerr(lun_err, err_out , 'read_mesh', &
            etb(fname),res)

    case('off')
       if ( this%nnodeincell .ne. 3) then
          if(res .ne. 0) rc = IOerr(lun_err, err_val , 'write_mesh', &
               ' format off works only for 2d and surface triangualtion')
       end if

       !
       ! write head
       !
       write(u_number,'(a)',iostat=res) 'OFF'
       if(res .ne. 0) rc = IOerr(lun_err, err_out , 'read_mesh', &
            etb(fname) // ' vars nnode',res)
       write(u_number,*,iostat=res) this%nnode, this%ncell, &
            this%nnode+this%ncell-2  
       if(res .ne. 0) rc = IOerr(lun_err, err_out , 'read_mesh', &
            etb(fname) // ' vars this%ncell, nnodeincell, mesh_type',res)

       !
       ! write coordinates
       !
       do j=1,this%nnode
          write(u_number,*,iostat=res) (this%coord(k,j),k=1,3)
          if(res .ne. 0) THEN
             write(rdwr,'(i5)') j
             str=etb(rdwr)//'/'
             rc = IOerr(lun_err, err_out , 'read_mesh', &
                  trim(etb(fname)) //&
                  ' type mesh member array coord at line '//&
                  trim(str),res)
          end if
       end do

       !
       ! write topol
       !
       do j=1,this%ncell
          write(u_number,*,iostat=res) this%nnodeincell, (this%topol(k,j)-1, k=1,this%nnodeincell)
          if(res .ne. 0) then
             write(rdwr,'(i5)') j
             str=etb(rdwr)//'/'
             write(rdwr,'(i5)') k
             str=trim(str)//etb(rdwr)
             rc = IOerr(lun_err, err_inp , 'write_mesh', &
                  trim(etb(fname)) // &
                  ' type mesh member array triang at line/col '//trim(str),res)
          end if
       end do


    end select


  end subroutine write_mesh

  !>-------------------------------------------------------------
  !> Build a graph describeing the cell to cell connection
  !> only those that share one edge/face
  !>
  !> M_{icell,inode}=1 if inode \in icell
  !>                 0 otherwise
  !> (procedure private for type mesh)
  !>
  !> usage:
  !>     call this%build_connection_matrix(lun_err, connection_matrix)
  !>
  !> where:
  !> \param[in] lun_err           -> integer. IO unit for error msg
  !> \param[in] connection_matrix -> type(spmatr). Sparse matrix with
  !>                                 M(icell,inode).ne. 0 if inode\in \icell
  !<-------------------------------------------------------------
  subroutine build_cellcell_graph(graph, lun_err, grid)
    use Globals
    use SparseMatrix
    implicit none
    !vars
    class(abs_simplex_mesh),        intent(inout) :: graph
    integer,                        intent(in   ) :: lun_err
    type(abs_simplex_mesh),target,  intent(in   ) :: grid
    ! local
    logical :: clean_memory
    type(abs_simplex_mesh), target :: local
    class(abs_simplex_mesh),pointer :: mesh

    select case (grid%nnodeincell)
    case (3)
       ! triangular cells

       ! create local copy if required
       if ( ( .not. allocated(grid%edge_iside) )  .or. ( .not. allocated(grid%bar_cell) ))  then
          local=grid
          call local%build_edge_connection(lun_err)
          call local%build_bar_cell()
          clean_memory = .true.
          mesh => local
       else
          clean_memory=.False.
          mesh => grid
       end if

       ! here we use that boundary edges are the first in edge_plist
       if (allocated(graph%topol)) call graph%kill(lun_err)
       call graph%init_from_data(lun_err, &
            mesh%ncell,mesh%nedge-mesh%nedge_bc, 2,'edge',&
            mesh%edge_plist(1:2,mesh%nedge_bc+1:mesh%nedge), mesh%bar_cell)
    case default
       write(lun_err,*) ' subroutine build_cellcell_graph works only for triangles'
    end select

    if (clean_memory) call mesh%kill(lun_err)
    
  end subroutine build_cellcell_graph

  
  !>-------------------------------------------------------------
  !> Build ncell*nnode sparse matrix M with 
  !>
  !> M_{icell,inode}=1 if inode \in icell
  !>                 0 otherwise
  !> (procedure private for type mesh)
  !>
  !> usage:
  !>     call this%build_connection_matrix(lun_err, connection_matrix)
  !>
  !> where:
  !> \param[in] lun_err           -> integer. IO unit for error msg
  !> \param[in] connection_matrix -> type(spmatr). Sparse matrix with
  !>                                 M(icell,inode).ne. 0 if inode\in \icell
  !<-------------------------------------------------------------
  subroutine build_connection_matrix(this, lun_err, connection_matrix)
    use Globals
    use SparseMatrix
    implicit none
    !vars
    class(abs_simplex_mesh), intent(in   ) :: this
    integer,                 intent(in   ) :: lun_err
    type(spmat),             intent(inout) :: connection_matrix
    ! local vars
    integer :: icell

    call connection_matrix%init(lun_err,&
         this%ncell, this%nnode, this%ncell* this%nnodeincell,&
         'csr', &
         is_symmetric = .False. )

    connection_matrix%ia(1)=1
    do icell = 1, this%ncell
       connection_matrix%ia(icell+1) = icell * this%nnodeincell + 1
       connection_matrix%ja(&
            connection_matrix%ia(icell):connection_matrix%ia(icell+1)-1 ) = &
            this%topol(1:this%nnodeincell,icell)
       call isort( this%nnodeincell,&
            connection_matrix%ja(&
            connection_matrix%ia(icell):connection_matrix%ia(icell+1)-1 ))
    end do
    connection_matrix%coeff = zero

  end subroutine build_connection_matrix

  !>-------------------------------------------------------------
  !> Build ncell*nnode sparse matrix M with 
  !>
  !> M_{icell,inode}=1 if inode \in icell
  !>                 0 otherwise
  !> (procedure private for type mesh)
  !>
  !> usage:
  !>     call this%build_nodenode_matrix(lun_err, connection_matrix)
  !>
  !> where:
  !> \param[in] lun_err           -> integer. IO unit for error msg
  !> \param[in] connection_matrix -> type(spmatr). Nnode x Nnode Sparse matrix with
  !>                                 M(i,j).ne.0 if node i is connected with 
  !>                                 node j. Self connection on the diagonal
  !>                                 can be included or not
  !<-------------------------------------------------------------
  subroutine build_nodenode_matrix(grid, lun_err,include_diagonal,adj_matrix)
    use Globals
    use SparseMatrix
    implicit none
    class(abs_simplex_mesh), target, intent(in) :: grid
    integer,      intent(in   ) :: lun_err
    logical,      intent(in   ) :: include_diagonal
    type(spmat),  intent(inout) :: adj_matrix
    !local
    logical :: rc
    logical :: clean_memory
    integer :: res
    integer :: iedge,inode,i1,i2,nnode,nterm,j
    integer, allocatable,dimension(:) :: iaw
    integer, allocatable,dimension(:) :: jaw
    integer, allocatable,dimension(:) :: count    
    class(abs_simplex_mesh),pointer :: this
    type(abs_simplex_mesh), target :: local

    if ( .not. allocated(grid%edge_iside))  then
       select type (grid)
       type is (abs_simplex_mesh)
          local = grid
       end select
       call local%build_edge_connection(lun_err)
       clean_memory = .true.
       this => local 
    else
       this => grid
    end if

    nnode = this%nnode
    if ( include_diagonal) then
       nterm = 2*this%nedge+this%nnode
    else
       nterm = 2*this%nedge
    end if

    iedge=nnode+1
    allocate(iaw(iedge),stat=res)
    if(res .ne. 0) rc = IOerr(lun_err, err_alloc, 'build_nodenode_matrix', &
         'work array iaw ',res)
    allocate(jaw(nterm),stat=res)
    if(res .ne. 0) rc = IOerr(lun_err, err_alloc, 'build_nodenode_matrix', &
         'work array jaw ',res)
    allocate(count(nnode),stat=res)
    if(res .ne. 0) rc = IOerr(lun_err, err_alloc, 'build_nodenode_matrix', &
         'work array  count',res)


    ! build ia pointer
    if ( include_diagonal) then
       count=1
    else
       count=0
    end if
    do iedge = 1,this%nedge
       i1=this%edge_iside(1,iedge)
       i2=this%edge_iside(2,iedge)
       count(i1) = count(i1)+1
       count(i2) = count(i2)+1
       iaw(this%edge_iside(1,iedge)) = iaw(this%edge_iside(1,iedge))+1 
       iaw(this%edge_iside(2,iedge)) = iaw(this%edge_iside(2,iedge))+1
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
       i1 = this%edge_iside(1,iedge)
       i2 = this%edge_iside(2,iedge)
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

    deallocate(iaw,jaw,count,stat=res)
    if(res .ne. 0) rc = IOerr(lun_err, err_dealloc, 'init_star', &
         'work array iaw jaw count',res)
    if(res .ne. 0) rc = IOerr(lun_err, err_dealloc, 'init_star', &
         'work array ia ja count',res)




    if (clean_memory) then
       call local%kill(lun_err)
    end if

  end subroutine build_nodenode_matrix

  !>-------------------------------------------------------------
  !> Build ncell*nnode sparse matrix M with 
  !>
  !> M_{icell,jcell}=1 if icell "touch" jcell
  !>                 0 otherwise
  !> (procedure private for type mesh)
  !>
  !> usage:
  !>     call this%build_nodenode_matrix(lun_err, connection_matrix)
  !>
  !> where:
  !> \param[in] lun_err           -> integer. IO unit for error msg
  !> \param[in] connection_matrix -> type(spmatr). Ncell x Ncell Sparse matrix with
  !>                                 M(i,j).ne.0 if cell i and cell j have one
  !>                                 edge face in common.
  !>                                 Self connection on the diagonal
  !>                                 can be included or not
  !<-------------------------------------------------------------
  subroutine build_cellcell_matrix(this, lun_err,include_diagonal,adj_matrix)
    use Globals
    use SparseMatrix
    implicit none
    class(abs_simplex_mesh), target, intent(in) :: this
    integer, intent(in   ) :: lun_err
    logical, intent(in   ) :: include_diagonal
    type(spmat), target, intent(inout) :: adj_matrix
    !local
    logical :: rc
    logical :: clean_memory
    integer :: res
    integer :: iedge,inode,icell,i1,i2
    integer :: nnode,nterm,ncell,nedge,nedge_bc
    integer, pointer :: ia(:)
    integer, pointer :: ja(:)
    integer, allocatable :: count(:)    
    class(abs_simplex_mesh),pointer :: grid
    type(abs_simplex_mesh), target  :: local

    !
    ! ensure to having edge_plist arrays
    !
    if (this%nnodeincell .ne. 3) then
       rc = IOerr(lun_err, err_val, 'select', &
         'curretn version supports only 2d-surface ')
    end if
    
    if ( allocated(this%edge_plist))  then
       grid => this
    else
       select type (this)
       type is (abs_simplex_mesh)
          local = this
       end select
       call local%build_edge_connection(lun_err)
       clean_memory = .true.
       grid => local 
    end if

    nnode    = grid%nnode
    ncell    = grid%ncell
    nedge    = grid%nedge
    nedge_bc = grid%nedge_bc

    if ( include_diagonal) then
       nterm = 2*(nedge-nedge_bc)+ncell
    else
       nterm = 2*(nedge-nedge_bc)
    end if

    !
    ! init matrix and temp arrays
    !
    call adj_matrix%init(lun_err,&
         ncell,ncell, nterm,&
         storage_system='csr',&
         is_symmetric=.true.)
    ia => adj_matrix%ia
    ja => adj_matrix%ja
    allocate(count(ncell),stat=res)
    if(res .ne. 0) rc = IOerr(lun_err, &
         err_alloc, 'build_nodenode_matrix', &
         'work array count',res)

    
    !
    ! build ia pointer
    !
    if ( include_diagonal) then
       count=1
    else
       count=0
    end if  
    do iedge = grid%nedge_bc+1,grid%nedge
       i1=grid%edge_plist(1,iedge)
       i2=grid%edge_plist(2,iedge)
       count(i1) = count(i1)+1
       count(i2) = count(i2)+1
    end do
    ia(1) = 1 
    do icell = 1, ncell
       ia(icell+1) =  ia(icell)+count(icell)
    end do

    !
    ! build ja
    !
    if ( include_diagonal ) then
       ! add self connection
       do icell = 1, grid%ncell
          ja(ia(icell)) = icell
       end do
       count = 1
    else
       count = 0
    end if   
    do iedge = grid%nedge_bc+1, grid%nedge
       i2 = grid%edge_plist(2,iedge)
       i1 = grid%edge_plist(1,iedge)
       ja( ia(i1)+count(i1) ) = i2 
       count( i1 ) = count( i1 ) + 1 
       ja( ia(i2)+count(i2) ) = i1
       count( i2 ) = count( i2 ) + 1 
    end do
    adj_matrix%is_sorted = .False.
    adj_matrix%coeff = one
   
    call adj_matrix%sort()

    !
    ! free memory
    !
    deallocate(count,stat=res)
    if(res .ne. 0) rc = IOerr(lun_err, err_dealloc, 'init_star', &
         'work array iaw jaw count',res)
    
    if (clean_memory) then
       call local%kill(lun_err)
    end if

  end subroutine build_cellcell_matrix


  !>-------------------------------------------------------------
  !> Part of the static constructor.
  !> (procedure public for type abs_simplex_mesh)
  !> Instantiate and initialize type abs_simplex_mesh
  !> from another via a selector of nodes or cells
  !>
  !> usage:
  !>     call this%read_mesh(lun_err, file2read)
  !>
  !> where:
  !> \param[in] lun_err     -> integer. IO unit for error msg
  !> \param[in] parent_mesh -> type(abs_mesh). Original mesh
  !> \param[in] ndata       -> integer. Lenght of selector
  !>                           Must be equal to ncel or nnode
  !> \param[in] selector    -> logical(dim =ndata) 
  !>                           If selector(i)=.False. remove node(i)/cell(i) 
  !<-------------------------------------------------------------
  subroutine init_selection(this, lun_err, parent_mesh, ndata, selector,&
       other_selector,&
       map_cells_new2old, map_nodes_new2old)
    use Globals
    class( abs_simplex_mesh), intent(inout) :: this
    integer,                  intent(in   ) :: lun_err
    class(abs_simplex_mesh),  intent(in   ) :: parent_mesh
    integer,                  intent(in   ) :: ndata
    logical,                  intent(in   ) :: selector(ndata)
    logical, optional,        intent(inout) :: other_selector(:)
    integer, optional,        intent(inout) :: map_cells_new2old(parent_mesh%ncell)
    integer, optional,        intent(inout) :: map_nodes_new2old(parent_mesh%nnode)
    
    ! local
    logical :: rc
    integer :: res
    integer :: inode, icell, i,j,iloc
    integer :: nnode, ncell,nnodeincell
    integer :: nnode_sub, ncell_sub
    integer, allocatable ::  topol_temp(:,:)
    integer, allocatable ::  index_inv(:)

    integer :: new_nnode,new_ncell
    integer, allocatable ::  new_node(:)
    integer, allocatable ::  inv_new_node(:)
    integer, allocatable ::  new_cell(:)
    integer, allocatable ::  inv_new_cell(:)
    logical, allocatable ::  selector_node(:)
    logical :: cell2preserve

    nnode       = parent_mesh%nnode
    ncell       = parent_mesh%ncell
    nnodeincell = parent_mesh%nnodeincell

    allocate(&
         topol_temp(nnodeincell+1,ncell),&
         index_inv(nnode),&
         stat=res)
    if(res .ne. 0) rc = IOerr(lun_err, err_alloc, 'select', &
         ' temp_array topol_temp,  index_inv',res)

    topol_temp = 0
    index_inv  = 0

    this%nnodeincell=parent_mesh%nnodeincell
    this%ambient_dimension=parent_mesh%ambient_dimension
    this%logical_dimension=parent_mesh%logical_dimension
    this%cell_type=parent_mesh%cell_type
    this%mesh_type=parent_mesh%mesh_type

    
    if ( ndata .eq. nnode) then
       !
       ! defied index_inv and new number of node
       !
       nnode_sub = 0
       do inode = 1, nnode
          if ( selector(inode) ) then
             nnode_sub               = nnode_sub + 1
             index_inv(inode)        = nnode_sub
          end if
       end do

       ncell_sub = 0
       do icell = 1, ncell
          !
          ! if any node in cell has been removed remove the cell
          !
          j=0
          cell2preserve = .True.
          do while ( (cell2preserve) .and. (j .lt. nnodeincell ) )
             j=j+1
             cell2preserve = cell2preserve .and. &
                  selector(parent_mesh%topol(j,icell))
          end do
          if ( cell2preserve) then
             ncell_sub = ncell_sub + 1
             topol_temp(1:nnodeincell,ncell_sub) = &
                  index_inv(parent_mesh%topol(1:nnodeincell,icell))
             ! store material
             topol_temp(1+nnodeincell,ncell_sub) = &
                  parent_mesh%topol(1+nnodeincell,icell)
          end if
          if (present(other_selector))  other_selector(icell) = cell2preserve
       end do

       if(nnode_sub .eq. 0) rc = IOerr(lun_err, err_val, 'selection_init', &
            ' null number of nodes, no output will be generated')
       if(ncell_sub .eq. 0) rc = IOerr(lun_err, err_val, 'selection_init', &
            ' null number of cell, no output will be generated')


       allocate( &
            this%coord(3,nnode_sub), &
            this%topol(nnodeincell+1,ncell_sub),&
            stat=res)



       !
       ! define new grid
       !  
       this%nnode=nnode_sub
       this%ncell=ncell_sub
       this%nnodeincell=parent_mesh%nnodeincell
       this%ambient_dimension=parent_mesh%ambient_dimension
       this%cell_type=parent_mesh%cell_type
       this%mesh_type=parent_mesh%mesh_type

       j=0
       do inode=1,nnode
          if ( selector(inode) ) then
             j=j+1
             this%coord(:,j)=parent_mesh%coord(:,inode)
          end if
       end do

       this%topol(:,1:ncell_sub) = topol_temp(:,1:ncell_sub) 

       !   
       ! fre memory
       !
       deallocate(&
            topol_temp,&
            index_inv,stat=res)

       if(res .ne. 0) rc = IOerr(lun_err, err_dealloc, 'init_selection', &
            ' temp_array topol_temp, index_inv',res)
    else if ( ndata .eq. ncell) then
       allocate(& 
            selector_node(nnode),&
            new_node(nnode),&
            inv_new_node(nnode),&
            new_cell(ncell),&
            inv_new_cell(ncell),&
            stat=res)
       if(res .ne. 0) rc = IOerr(lun_err, err_alloc, 'init_selection', &
            ' temp_array ',res)


       selector_node = .False.
       new_ncell= 0
       new_cell(:) = 0
       inv_new_cell(:) = 0
       do icell = 1, parent_mesh%ncell
          if (selector(icell) ) then
             new_ncell = new_ncell + 1
             new_cell(icell)  = new_ncell 
             inv_new_cell(new_ncell) = icell
          end if
          do iloc= 1, parent_mesh%nnodeincell
             inode = parent_mesh%topol(iloc,icell) 
             selector_node (inode) = selector_node (inode) .or. selector(icell)
          end do
       end do

       new_nnode = 0
       new_node(:) = 0
       inv_new_node(:) = 0
       do inode = 1, parent_mesh%nnode
          if ( selector_node (inode) ) then
             new_nnode=new_nnode+1
             new_node(inode) = new_nnode
             inv_new_node(new_nnode) = inode
          end if
       end do

       this%nnode = new_nnode
       this%ncell = new_ncell
      
       allocate( &
            this%coord(3,new_nnode),&
            this%topol(this%nnodeincell+1,new_ncell),&
            stat=res)

       do inode = 1, new_nnode
          this%coord(:,inode) = parent_mesh%coord(:, inv_new_node(inode))
       end do

       do icell = 1, new_ncell
          ! copy topol in new nodes ordering
          do iloc = 1, this%nnodeincell
             this%topol(iloc,icell) = &
                  new_node(parent_mesh%topol(iloc, inv_new_cell(icell)))
          end do
          ! copy material
          iloc = this%nnodeincell +1 
          this%topol(iloc,icell) = parent_mesh%topol(iloc, inv_new_cell(icell) )
       end do

       if ( present (other_selector) ) other_selector = selector_node

       if ( present ( map_cells_new2old) ) map_cells_new2old = inv_new_cell      
       if ( present ( map_nodes_new2old) ) map_nodes_new2old = inv_new_node
       
       deallocate(& 
            selector_node,&
            new_node,&
            inv_new_node,&
            new_cell,&
            inv_new_cell,&
            stat=res)
       if(res .ne. 0) rc = IOerr(lun_err, err_dealloc, 'init_selection', &
            ' temp_array ',res)



    else
       rc = IOerr(lun_err, err_val, 'init_selection', &
            ' selector not compatible with mesh ndata =', ndata)
    end if


  end subroutine init_selection

  subroutine mesh_constructor_from_cell_selection(this, lun_err, parent_mesh, selector,&
       other_selector,&
       map_cells_new2old,map_nodes_new2old)
    use Globals
    type(abs_simplex_mesh), intent(inout) :: this
    integer,                intent(in   ) :: lun_err
    type(abs_simplex_mesh), intent(in   ) :: parent_mesh
    logical,                intent(in   ) :: selector(parent_mesh%ncell)
    logical,         intent(inout) :: other_selector(parent_mesh%nnode)
    integer,         intent(inout) :: map_cells_new2old(parent_mesh%ncell)
    integer,         intent(inout) :: map_nodes_new2old(parent_mesh%nnode)

    call this%init_selection(lun_err, parent_mesh, this%ncell, selector,&
         other_selector,&
         map_cells_new2old,map_nodes_new2old)

  end subroutine mesh_constructor_from_cell_selection


  !>------------------------------------------------------------------
  !> procedure for projection over a rifened mesh vreal var. data
  !> (public procedure for type mesh)
  !>
  !> usage:    call  proj_subgrid(subgrid,data,data_avg)
  !>
  !> where:
  !> \param[in ] data       -> real(ncell_parent) data to be averaged
  !> \param[out] data_prj   -> real(ncell_sub)    data averaged
  !>------------------------------------------------------------------
  subroutine proj_subgrid(subgrid, data, data_prj)
    use Globals
    implicit none
    class(abs_simplex_mesh),       intent(in   ) :: subgrid
    real(kind=double), intent(in   ) :: data(subgrid%ncell_parent)
    real(kind=double), intent(inout) :: data_prj(subgrid%ncell)
    !local 
    integer :: icell_sub,i

    do icell_sub = 1, subgrid%ncell
       i=subgrid%cell_parent(icell_sub)
       data_prj(icell_sub) = data(i)
       !write(*,*) icell_sub,i,data(i), data_prj(icell_sub)
    end do

  end subroutine proj_subgrid

  !>------------------------------------------------------------------
  !> procedure for projection over a rifened mesh vreal var. data
  !> (public procedure for type mesh)
  !>
  !> usage:    call  projnode_subgrd(subgrid,data,data_avg)
  !>
  !> where:
  !> \param[in ] data       -> real(nnode_parent) data to be averaged
  !> \param[out] data_prj   -> real(nnode_sub) data averaged
  !>------------------------------------------------------------------
  subroutine projnode_subgrid(subgrid, data, data_prj)
    use Globals
    implicit none
    class(abs_simplex_mesh),       intent(in   ) :: subgrid
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
  !> procedure for projection over a rifened mesh vreal var. data
  !> (public procedure for type mesh)
  !>
  !> usage:    call  proj_subgrid(subgrid,data,data_avg)
  !>
  !> where:
  !> \param[in ] data       -> real(ncell_parent) data to be averaged
  !> \param[out] data_prj   -> real(ncell_sub)    data averaged
  !>------------------------------------------------------------------
  subroutine avg_cell_subgrid(subgrid, grid, data, data_prj)
    use Globals
    implicit none
    class(abs_simplex_mesh),  intent(in   ) :: subgrid
    class(abs_simplex_mesh),  intent(in   ) :: grid
    real(kind=double),        intent(in   ) :: data(subgrid%ncell)
    real(kind=double),        intent(inout) :: data_prj(grid%ncell)
    !local 
    integer :: icell,icell_sub

    data_prj = zero

    do icell_sub = 1, subgrid%ncell 
       icell = subgrid%cell_parent(icell_sub)
       data_prj(icell) = data_prj(icell) + &
            data(icell_sub) * subgrid%size_cell(icell_sub)  / grid%size_cell(icell) 
    end do

  end subroutine avg_cell_subgrid


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
    class(abs_simplex_mesh),    intent(inout) :: this
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
    this%grid_level = parent_level

    ! read  array dimension and check that the mesh is consistent
    read(u_number,*,iostat=res) nnode_parent_read,  subnnode_read 
    if(res .ne. 0) rc = IOerr(lun_err, err_inp , 'read_grid_subgrid_relation', &
         etb(fname) // ' vars nnode_parent_read,  subnnode_read ',res)
    if ( this%nnode >0) then
       if ( subnnode_read .ne. this%nnode) then
          write(lun_err,*) ' Mismatch node number'
          write(lun_err,*) ' subgrid nnode = ', this%nnode
          write(lun_err,*) ' subnnode read  = ', subnnode_read 
          stop
       end if
    else
       this%nnode = subnnode_read 
    end if

    read(u_number,*,iostat=res) ncell_parent_read,  subncell_read 
    if(res .ne. 0) rc = IOerr(lun_err, err_inp , 'read_grid_subgrid_relation', &
         etb(fname) // ' vars ncell_parent_read, subncell_read',res)
    if ( this%ncell >0) then
       if ( subncell_read .ne. this%ncell) then
          write(lun_err,*) ' Mismatch node number'
          write(lun_err,*) ' subgrid ncell = ', this%ncell
          write(lun_err,*) ' subncell read = ', subncell_read 
          stop
       end if
    else
       this%ncell = subncell_read
       select case( this%ncell/ncell_parent_read)
       case (4)
          this%logical_dimension=2
       case (8)
          this%logical_dimension=2
       case (2)
          this%logical_dimension=1
       end select
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
    class(abs_simplex_mesh),    intent(in) :: this
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
    write(u_number,*,iostat=res) this%grid_level, ' ! grid level=number of refinements between two grids '
    if(res .ne. 0) rc = IOerr(lun_err, err_out , 'write_grid_subgrid_relation', &
         etb(fname) // ' vars grid_level',res)

    if (this%grid_level.eq.0) then
       rc=IOerr(lun_err, wrn_out , &
            'write_grid_subgrid_relation', &
            etb(fname) // ' grid_level 0: no parent info',res)
       return
    end if

    ! Write array dimensions
    write(u_number,*,iostat=res) this%nnode_parent, this%nnode, ' ! number of node coaser grid,  number of nnode finer grid '
    if(res .ne. 0) rc = IOerr(lun_err, err_out , 'write_grid_subgrid_relation', &
         etb(fname) // ' vars nnode',res)

    write(u_number,*,iostat=res) this%ncell_parent, this%ncell, ' ! number of cell coaser grid,  number of ncell finer grid '
    if(res .ne. 0) rc = IOerr(lun_err, err_out, 'write_grid_subgrid_relation', &
         etb(fname) // ' vars ncell',res)

    ! Write node parent
    write(u_number,*,iostat=res) (this%node_parent(k,1),k=1,2), ' ! n1, n2 : first and second nodes parents '
    if(res .ne. 0) then
       write(rdwr,'(i5)') 1
       str=trim(adjustl(rdwr))//'/'
       rc = IOerr(lun_err, err_out , 'write_grid_subgrid_relation', &
            trim(etb(fname)) //&
            ' type mesh member array coord at line '//&
            trim(str),res)
    end if
    do inode=2,this%nnode
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

    !
    ! Write cell parent
    ! 
    write(u_number,*,iostat=res) this%cell_parent(1), ' ! index of cell in coaser grid '
    if(res .ne. 0) then
       write(rdwr,'(i6)') 1
       str=trim(adjustl(rdwr))
       rc = IOerr(lun_err, err_out , 'write_grid_subgrid_relation', &
            trim(etb(fname)) // &
            ' type mesh member array tetra at tetra nmb: '//trim(str),res)
    end if
    do icell=2,this%ncell
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



  !>------------------------------------------------------------------
  !> Procedure for costrustion of size_cell member for graph, 
  !> triangular and tetrahedral grid.
  !> (public procedure for type mesh)
  !> It allocates if necessary.
  !>
  !> usage:    call  var%build_size_cell(subgrid,data,data_avg)
  !>
  !> where:
  !> \param[in] lun_err  -> integer. unit for err msg.
  !>------------------------------------------------------------------
  subroutine build_size_cell(this,lun_err)
    use Globals
    implicit none
    class(abs_simplex_mesh),  intent(inout) :: this
    integer,                  intent(in   ) :: lun_err

    !local
    logical :: rc
    integer :: res
    integer :: icell
    real(kind=double) :: dnrm2
    real(kind=double) :: v1(3),v2(3),v3(3)

    if (  .not. allocated(this%size_cell) ) then
       allocate( this%size_cell(this%ncell) ,stat=res) 
       if (res .ne. 0) &
            rc = IOerr(lun_err, err_inp , ' build size cell ', &
            ' type abs_simplex_mesh member size_cell ', res)
    end if

    select case ( this%nnodeincell) 
    case (4)
       select case ( this%cell_type) 
       case ('tetrahedron') 
          do icell = 1,this%ncell
             v1 = this%coord(1:3,this%topol(1,icell))-this%coord(1:3,this%topol(4,icell))
             v2 = this%coord(1:3,this%topol(2,icell))-this%coord(1:3,this%topol(4,icell))
             v3 = this%coord(1:3,this%topol(3,icell))-this%coord(1:3,this%topol(4,icell))
             this%size_cell(icell) = volume(v1,v2,v3)
          end do
       case ('square')
          !
          ! assuming to number 
          !  4   3
          !
          !  1   2 
          ! 
          do icell = 1,this%ncell
             v1(1) = this%coord(1,this%topol(2,icell))-this%coord(1,this%topol(1,icell))
             v2(1) = this%coord(2,this%topol(4,icell))-this%coord(2,this%topol(1,icell))
             this%size_cell(icell) = abs(v1(1)) * abs(v2(1))
          end do
       case default
          rc = IOerr(lun_err, err_inp , ' build size cell ', &
               ' cell type'//etb(this%cell_type)//' not supported ')
       end select


    case (3)
       do icell = 1,this%ncell
          v1 = this%coord(1:3,this%topol(1,icell))-this%coord(1:3,this%topol(3,icell))
          v2 = this%coord(1:3,this%topol(2,icell))-this%coord(1:3,this%topol(3,icell))
          this%size_cell(icell) = area(v1,v2)
       end do
    case (2)
       do icell = 1,this%ncell
          v1 = this%coord(1:3,this%topol(1,icell))-this%coord(1:3,this%topol(2,icell))
          this%size_cell(icell) = dnrm2(2,v1,1)
       end do

    end select

  contains
    function volume(v1,v2,v3) result(vol)
      use Globals
      implicit none
      real(kind=double), intent(in) :: v1(3),v2(3),v3(3)
      real(kind=double) :: vol
      ! local
      real(kind=double) :: v2_cross_v3(3)
      real(kind=double) :: ddot

      v2_cross_v3=cross(v2,v3)
      vol=abs(ddot(3,v1,1,v2_cross_v3,1))*onesixth

    end function volume

    function area(v1,v2) result(ar)
      use Globals
      implicit none
      real(kind=double), intent(in) :: v1(3),v2(3)
      real(kind=double) :: ar
      ! local
      real(kind=double) :: v1_cross_v2(3)
      real(kind=double) :: dnrm2



      v1_cross_v2 = cross(v1,v2)
      ar=dnrm2(3,v1_cross_v2,1)*onehalf

    end function area


  end subroutine build_size_cell

  !>------------------------------------------------------------------
  !> Procedure for costrustion of normal_cell member for  
  !> triangular grid or surface grid.
  !> (public procedure for type mesh)
  !> It allocates if necessary.
  !>
  !> usage:    call  var%build_normal_cell(subgrid,data,data_avg)
  !>
  !> where:
  !> \param[in] lun_err  -> integer. unit for err msg.
  !>------------------------------------------------------------------
  subroutine build_bar_cell(this,lun_err)
    use Globals
    implicit none
    class(abs_simplex_mesh), intent(inout) :: this
    integer, optional,       intent(in   ) :: lun_err
    !local
    logical :: rc
    integer :: res,lun
    integer :: icell,iloc
    
    if (present(lun_err)) then
       lun=lun_err
    else
       lun=0
    end if

    if (  .not. allocated(this%bar_cell) ) then
       allocate( this%bar_cell(3,this%ncell) ,stat=res) 
       if (res .ne. 0) &
            rc = IOerr(lun_err, err_inp , ' build size cell ', &
            ' type abs_simplex_mesh member bar_cell ', res)
    end if
    !write(*,*)  'inside',allocated(this%bar_cell)

    this%bar_cell=zero
    do icell=1,this%ncell
       do iloc=1,this%nnodeincell
          this%bar_cell(:,icell)=this%bar_cell(:,icell)+&
               this%coord(:,this%topol(iloc,icell))
       end do
       this%bar_cell(:,icell)=this%bar_cell(:,icell)*one/this%nnodeincell
    end do
       
    

  end subroutine build_bar_cell

  !>------------------------------------------------------------------
  !> Procedure to compute nodal area/volume for graph, 
  !> triangular and tetrahedral grid.
  !> (public procedure for type mesh)
  !> It allocates if necessary.
  !>
  !> usage:    call  var%build_size_node(subgrid,data,data_avg)
  !>
  !> where:
  !> \param[in] lun_err  -> integer. unit for err msg.
  !>------------------------------------------------------------------
  subroutine build_size_node(this,lun_err)
    use Globals
    implicit none
    class(abs_simplex_mesh),  intent(inout) :: this
    integer,                  intent(in   ) :: lun_err
    ! local vars
    integer :: icell, inode, iloc
    real(kind=double) :: loc_val

    this%size_node = zero
    do icell=1,this%ncell
       loc_val = this%size_cell(icell)/this%nnodeincell

       do iloc=1,this%nnodeincell
          inode = this%topol(iloc,icell)

          this%size_node(inode) =this%size_node(inode) + loc_val
       end do
    end do

  end subroutine build_size_node


  !>------------------------------------------------------------------
  !> Procedure for costrustion of normal_cell member for  
  !> triangular grid or surface grid.
  !> (public procedure for type mesh)
  !> It allocates if necessary.
  !>
  !> usage:    call  var%build_normal_cell(subgrid,data,data_avg)
  !>
  !> where:
  !> \param[in] lun_err  -> integer. unit for err msg.
  !>------------------------------------------------------------------
  subroutine build_normal_cell(this,lun_err)
    use Globals
    implicit none
    class(abs_simplex_mesh),  intent(inout) :: this
    integer,                  intent(in   ) :: lun_err

    !local
    logical :: rc
    integer :: res
    integer :: icell
    integer :: n1,n2,n3
    real(kind=double) :: vec12(3),vec13(3),cross_prod(3)
    real(kind=double) :: dnrm2



    if ( this%nnodeincell .eq. 3) then
       if (  .not. allocated(this%normal_cell) ) then
          allocate( this%normal_cell(3,this%ncell) ,stat=res) 
          if (res .ne. 0) &
               rc = IOerr(lun_err, err_inp , ' build size cell ', &
               ' type abs_simplex_mesh member size_cell ', res)
       end if
       do icell=1,this%ncell
          n1 = this%topol(1,icell)
          n2 = this%topol(2,icell)
          n3 = this%topol(3,icell)
          vec12 = this%coord(:,n2)-this%coord(:,n1)
          vec13 = this%coord(:,n3)-this%coord(:,n1)
          cross_prod = cross(vec12,vec13)
          this%normal_cell(:,icell) = cross_prod/dnrm2(3,cross_prod,1)
       end do
    end if

  end subroutine build_normal_cell



  !>-------------------------------------------------------------
  !> Free memory procedure, deallocate all allocated members
  !> (procedure public for type abs_simplex_mesh)
  !>
  !> usage:
  !>     call 'varr'%info(lun)
  !>
  !> where:
  !> \param[in] lun  -> integer. IO unit 
  !<-------------------------------------------------------------
  subroutine kill_mesh(this, lun)
    use Globals
    implicit none
    !vars
    class(abs_simplex_mesh), intent(inout) :: this
    integer,                 intent(in   ) :: lun
    !local
    logical :: rc
    integer :: res
    type(abs_simplex_mesh) :: clean_mesh

    !
    ! we deallocate and reset all variables
    ! using '=' operator. 
    !
    select type (this)
    type is (abs_simplex_mesh)
       this=clean_mesh
    end select

    if ( allocated(this%size_cell) ) then
       deallocate( this%size_cell ,stat=res) 
       if (res .ne. 0) &
            rc = IOerr(lun, err_dealloc , 'kill_mesh ', &
            ' type abs_simplex_mesh member size_cell ', res)
    end if

    if ( allocated(this%normal_cell) ) then
       deallocate( this%normal_cell ,stat=res) 
       if (res .ne. 0) &
            rc = IOerr(lun, err_dealloc , 'kill_mesh ', &
            ' type abs_simplex_mesh member normal_cell ', res)
    end if

    if ( allocated(this%node_parent) ) then
       deallocate( this%node_parent ,stat=res) 
       if (res .ne. 0) &
            rc = IOerr(lun, err_dealloc , 'kill_mesh ', &
            ' type abs_simplex_mesh member node_parent ', res)
    end if


    if ( allocated(this%cell_parent) ) then
       deallocate( this%cell_parent ,stat=res) 
       if (res .ne. 0) &
            rc = IOerr(lun, err_dealloc , 'kill_mesh ', &
            ' type abs_simplex_mesh member cell_parent ', res)
    end if


    if ( allocated(this%coord) ) then
       deallocate( this%coord ,stat=res) 
       if (res .ne. 0) &
            rc = IOerr(lun, err_dealloc , 'kill_mesh ', &
            ' type abs_simplex_mesh member coord ', res)
    end if


    if ( allocated(this%topol) ) then
       deallocate( this%topol ,stat=res) 
       if (res .ne. 0) &
            rc = IOerr(lun, err_dealloc , 'kill_mesh ', &
            ' type abs_simplex_mesh member topol ', res)
    end if

     if ( allocated(this%edge_iside) ) then
       deallocate( this%edge_iside ,stat=res) 
       if (res .ne. 0) &
            rc = IOerr(lun, err_dealloc , 'kill_mesh ', &
            ' type abs_simplex_mesh iside mesh ', res)
    end if


    this%nnode = 0
    this%ncell = 0
    this%nedge = 0

    this%is_initialized=.False.
  end subroutine kill_mesh

  subroutine mesh_destructor(this, lun)
    use Globals
    implicit none
    !vars
    type(abs_simplex_mesh), intent(inout):: this
    integer,                intent(in   ) :: lun

    call this%kill(lun)
  end subroutine mesh_destructor

  !>------------------------------------------------------------------
  !> Evaluates the p-norm of functions piecewise constant on triangles
  !>------------------------------------------------------------------
  function normp_cell(this,power,data) result(total)
    use Globals
    implicit none
    class(abs_simplex_mesh), intent(in) :: this
    real(kind=double),       intent(in) :: power
    real(kind=double),       intent(in) :: data(this%ncell)
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
  subroutine build_face_connection(this, lun)
    use Globals
    implicit none
    !  External variables 
    class(abs_simplex_mesh), intent(inout) :: this
    integer,     intent(in) :: lun
    ! Local variables
    ! work arrays
    integer, allocatable :: topol(:,:)
    integer, allocatable :: faces(:,:), point(:)
    integer, allocatable :: neigh(:,:),face_cnc(:,:)
    integer, allocatable :: iside(:,:),face_plist(:,:)
    integer, allocatable :: perm_faces(:),perminv_faces(:)
    !  scalars and static arrays
    integer :: perm(3,4)
    data perm/1,2,3, 1,2,4, 1,3,4, 2,3,4/
    integer :: nface,nface_bc,ncell
    integer :: res
    integer :: icell,jtetra,ktetra,jj,kk,iloc
    integer :: iface,jface,kface
    integer :: i,j, indx,isgn, itemp
    integer :: temp(4)
    logical :: rc

	 if (this%logical_dimension.ne.3) then
		 write(*,*) 'ERROR in procedure build_face_connection: &
			 can be used for tetrahedral element'
		 STOP
	 end if

    ! Assign local scalar variabls
    ncell = this%ncell
    !
    ! Allocate local and work arrays
    ! 
    allocate(topol(4,ncell),faces(5,4*ncell),point(4*ncell),face_cnc(4,ncell),&
		 face_plist(2,4*ncell),iside(3,4*ncell),neigh(4,ncell),stat=res)
    if(res .ne. 0) rc = IOerr(lun, err_alloc, 'connection_faces', &
         'work and local arrays topol, faces, point, face_cnc, face_plist, iside, neigh',res)
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
    ! eliminate double faces and build face_cnc
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
    !  build the neigh vector using face_plist as a work array
    !
    do icell = 1,nface
       point(icell) = 1
    end do

      ! set face_cnc and neigh to 0
      face_cnc = 0
      neigh = 0

      do icell = 1,nface
         iside(1,icell) = faces(1,icell)
         iside(2,icell) = faces(2,icell)
         iside(3,icell) = faces(3,icell)
         face_plist(1,icell) = faces(4,icell)
         face_plist(2,icell) = faces(5,icell)
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
            if (face_cnc(kk,jtetra) .eq. 0) then
               face_cnc(kk,jtetra) = icell
               kk = 5
            end if
            kk = kk + 1
         end do
         if (ktetra .gt. 0) then
            kk = 1
            do while (kk .le. 4)
               if (face_cnc(kk,ktetra) .eq. 0) then
                  face_cnc(kk,ktetra) = icell
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
         if(face_plist(1,iface)*face_plist(2,iface).eq.0) then
            ! for the faces at the boundary
            kface=kface+1
            perm_faces(kface)=iface
         end if
      end do
      nface_bc=kface
      do iface=1,nface
         if(face_plist(1,iface)*face_plist(2,iface).ne.0) then
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
      if ( .not. allocated(this%face_cnc) ) then
         allocate(this%face_cnc(4,ncell),stat=res)
         if(res .ne. 0) rc = IOerr(lun, err_alloc, 'connection_faces', &
              ' type mesh member array face_cnc',res)
      end if

      if ( .not. allocated(this%face_iside) ) then
         allocate(this%face_iside(3,nface),stat=res)
         if(res .ne. 0) rc = IOerr(lun, err_alloc, 'connection_faces', &
              ' type mesh member array face_iside',res)
      end if
      if ( .not. allocated(this%face_plist) ) then
         allocate(this%face_plist(2,nface),stat=res)
         if(res .ne. 0) rc = IOerr(lun, err_alloc, 'connection_faces', &
              ' type mesh member array face_plist',res)
      end if

      this%nface = nface
      this%nface_bc = nface_bc

      do iface=1,nface
         this%face_iside(:,iface)=iside(:,perm_faces(iface))
         this%face_plist(:,iface)=face_plist(:,perm_faces(iface))
      end do
      do icell=1,ncell
         do iloc=1,4
            this%face_cnc(iloc,icell)=perminv_faces(face_cnc(iloc,icell))
         end do
      end do
      
      this%neigh_face=neigh    
      
    ! deallocate local and work arrays
    deallocate(faces,point,topol,face_cnc,face_plist,iside,neigh,perm_faces,perminv_faces,stat=res)
    if(res .ne. 0) rc = IOerr(lun, err_dealloc, 'connection_faces', &
    'work arrays faces, point, topol, face_cnc, face_plist, iside, neigh, perm_faces, perminv_faces, perm',res)

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

  end subroutine build_face_connection


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
  subroutine build_edge_connection(this, lun)
    use Globals
    implicit none
    !  External variables 
    class(abs_simplex_mesh), intent(inout) :: this
    integer,     intent(in   ) :: lun
    ! Local variables
    ! work arrays
    integer, allocatable :: topol(:,:)
    integer, allocatable :: edges(:,:)
    integer, allocatable :: point(:),point2(:),perm_edges(:),perminv_edges(:)
    integer, allocatable :: edge_iside(:,:),edge_cnc(:,:),edge_plist(:,:)
    !  scalars and static arrays
    integer :: perm(2,6)
    !data perm/1,2, 1,3, 1,4, 2,3, 2,4, 3,4/
    integer :: ncell, nedge_temp, count_edge
    integer :: res
    integer :: icell,jcell,kcell,iedge,jedge,kedge,ii,jj,kk,iloc
    integer :: i,j, indx,isgn, itemp
    integer :: temp(3)
    logical :: rc
    integer :: nedgeincell,nnodeincell
    integer , allocatable :: permutation(:,:)

    ! Assign local scalar variables
    ncell=this%ncell
    nnodeincell = this%nnodeincell
    if ( this%nnodeincell .eq. 4) then
       nedgeincell = 6
    else if (this%nnodeincell .eq. 3) then
       nedgeincell = 3
    end if
    allocate (permutation(2,nedgeincell))
    if ( this%nnodeincell .eq. 4) then
       permutation(:,1) = (/1,2/)
       permutation(:,2) = (/1,3/)
       permutation(:,3) = (/1,4/)
       permutation(:,4) = (/2,3/)
       permutation(:,5) = (/2,4/) 
       permutation(:,6) = (/3,4/)
    else if (this%nnodeincell .eq. 3) then
       permutation(:,1) = (/1,2/)
       permutation(:,2) = (/2,3/)
       permutation(:,3) = (/1,3/)
    end if

    
    nedge_temp=nedgeincell*ncell
        
    !
    ! Allocate local and work arrays
    !
    allocate(&
         topol(nnodeincell+1,ncell),&
         edges(4,nedge_temp),&
         stat=res)
    if(res .ne. 0) rc = IOerr(lun, err_alloc, 'connection_edges', &
         'work and local arrays topol, edges ',res)
    !
    ! Generate edges vector
    !
    topol=this%topol

    do icell = 1,ncell
       ! sort the element nodes in increasing order
       call isort(nnodeincell,topol(1:nnodeincell,icell))
       do jcell = 1,nedgeincell
          kk = nedgeincell*(icell-1) + jcell
          edges(1,kk) = topol(permutation(1,jcell),icell)
          edges(2,kk) = topol(permutation(2,jcell),icell)
          edges(3,kk) = icell
          edges(4,kk) = 0
       enddo
    enddo
    !
    ! sort in lexicographic order the edges vector
    !

    !  Initialize.
    i = 0
    indx = 0
    isgn = 0
    j = 0
    do 
       call global_heapsort(nedge_temp, indx, i,j,isgn)
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



    !
    ! eliminate double faces and build side_cnc
    !
    count_edge = 1
    edges(4,1) = count_edge
    do iedge=2,nedge_temp
       if ( (edges(1,iedge) .eq. edges(1,iedge-1)) .and. &
            (edges(2,iedge) .eq. edges(2,iedge-1))) then
          edges(4,iedge) = count_edge
       else
          count_edge = count_edge + 1
          edges(4,iedge) = count_edge
       end if
    end do

    !
    ! Store real number of edges
    !
    this%nedge=count_edge 

    !
    ! build iside_cnc
    !
    if ( .not. allocated(this%edge_iside) ) then
       allocate(this%edge_iside(2,this%nedge),stat=res)
       if(res .ne. 0) rc = IOerr(lun, err_alloc, &
            'connection_edges', &
            ' type mesh member array edge_iside',res)
    end if
    
    !
    ! build iside
    !
    this%edge_iside(1:2,1) = edges(1:2,1)
    count_edge = 1
    do iedge=2,nedge_temp
       if (edges(4,iedge) .ne. count_edge) then
          count_edge = count_edge + 1
          this%edge_iside(:,count_edge) = edges(1:2,iedge)
       end if
    end do

    !
    ! build side connection
    !
    if ( .not. allocated(this%edge_cnc) ) then       
       allocate(this%edge_cnc(nedgeincell,this%ncell),stat=res)
       if(res .ne. 0) rc = IOerr(lun, err_alloc, &
            ' connection_edges ', &
            ' type mesh member array edge_cnc',res)
    end if
    
    !
    ! build edge_cnc
    ! initialize edge_cnc to zero
    !
    this%edge_cnc = 0
    do iedge=1,nedge_temp
       icell = edges(3,iedge)
       kk = 1
       do while (kk .le. nedgeincell)
          if (this%edge_cnc(kk,icell) .eq. 0) then
             this%edge_cnc(kk,icell) = edges(4,iedge)
             kk = nedgeincell+1
          end if
          kk = kk + 1
       end do
    end do

    !
    ! for 2d or surface grid build 
    ! edge_plist  : list of left-right triangles 
    ! neight : list of surrundong triangles for each triangles
    ! Moreover, permute edges thus that boundary edges are first
    !
    if (nnodeincell .eq. 3) then
       !
       !  build the neigh vector 
       !
       allocate(&
            point(this%ncell),&
            point2(this%nedge),&
            perm_edges(this%nedge),&
            perminv_edges(this%nedge),&
            stat=res)
       if(res .ne. 0) rc = IOerr(lun, err_alloc, 'connection_edges', &
            ' temporary array point point2 perm_edges, perminv_edges ', res)

       if ( .not. allocated(this%edge_plist) ) then
          allocate(this%edge_plist(2,this%nedge),stat=res)
          if(res .ne. 0) rc = IOerr(lun, err_alloc,&
               'connection_edges', &
               ' type mesh member edge_plist',res)
       end if

       if ( .not. allocated(this%neigh_cell) ) then
          allocate(this%neigh_cell(nedgeincell,this%ncell),stat=res)
          if(res .ne. 0) rc = IOerr(lun, err_alloc, &
               'connection_edges', &
               ' type mesh member neigh_cell',res)
       end if

       point=1
       point2=1
       this%neigh_cell=0
       this%edge_plist=0
       do i = 1, nedge_temp
          icell = edges(3,i)
          iedge = edges(4,i)

          !
          ! all edge
          !
          this%neigh_cell(point(icell),icell) = iedge
          point(icell) = point(icell) + 1
          
          !
          ! add cell
          !
          this%edge_plist(point2(iedge),iedge) = icell
          point2(iedge) = point2(iedge) + 1
       end do


       !
       ! find boundary edges
       !
       kedge=0
       do iedge=1,this%nedge
          if(this%edge_plist(1,iedge)*this%edge_plist(2,iedge).eq.0) then
             ! for the edges at the boundary
             kedge=kedge+1
             perm_edges(kedge)=iedge
          end if
       end do
       this%nedge_bc=kedge
       do iedge=1,this%nedge
          if(this%edge_plist(1,iedge)*this%edge_plist(2,iedge).ne.0) then
             ! for the internal edges
             kedge=kedge+1
             perm_edges(kedge)=iedge
          end if
       end do
       do iedge=1,this%nedge
          perminv_edges(perm_edges(iedge))=iedge
       end do



       !
       ! sort all edge quanties
       !
       allocate(&
            edge_iside(2,this%nedge),&
            edge_plist(2,this%nedge),&
            edge_cnc(3,this%nedge),&
            stat=res)
       if(res .ne. 0) rc = IOerr(lun, err_alloc, 'connection_edges', &
            ' temporary array point perm_edges, perminv_edges ', res)

       edge_iside= this%edge_iside
       edge_plist     = this%edge_plist
       edge_cnc  = this%edge_cnc

       
       do iedge=1,this%nedge
          this%edge_iside(:,iedge)=edge_iside(:,perm_edges(iedge))
          this%edge_plist(:,iedge)=edge_plist(:,perm_edges(iedge))
       end do
       do icell=1,ncell
          do iloc=1,nedgeincell
             this%edge_cnc(iloc,icell)=perminv_edges(edge_cnc(iloc,icell))
          end do
       end do


       !
       ! free memory
       !
       deallocate(&
            edge_iside,&
            edge_plist,&
            edge_cnc,&
            stat=res)
       if(res .ne. 0) rc = IOerr(lun, err_dealloc, 'connection_edges', &
            ' temporary array edge_iside edge_plist edge_cnc ', res)
       
       deallocate(point,point2,perm_edges,perminv_edges,stat=res)
       if(res .ne. 0) rc = IOerr(lun, err_dealloc, 'connection_edges', &
            ' temporary array point perm_edges, perminv_edges ', res)

    end if

  end subroutine build_edge_connection

  !>-------------------------------------------------------------
  !> Build list of boundary nodes
  !> (private procedure for type mesh )
  !> 
  !> usage:
  !>     call 'var'%build_nodebc(lun)
  !>
  !> where:
  !> \param[in], optional, lun_err -> integer. I/O unit for error message output
  !<-------------------------------------------------------------
  subroutine build_nodebc(this,lun_err)
    use KindDeclaration
    use Globals
    implicit none
    class(abs_simplex_mesh), intent(inout) :: this
    integer, optional,       intent(in   ) :: lun_err
    !local
    logical :: rc
    integer :: res,lun
    integer :: i,j,m
    integer, allocatable :: temp(:), temp1(:), temp2(:)


    if (present(lun_err)) then
       lun=lun_err
    else
       lun=0
    end if

    write(*,*) this%logical_dimension
    if (this%logical_dimension.eq.2) then
		 if (allocated(this%edge_iside)) then
			 this%nnode_bc=this%nedge_bc
			 allocate(temp(2*this%nedge_bc),stat=res)
			 if (res.ne.0)  rc = IOerr(lun, err_alloc, 'build_nodebc', &
					' temporary array temp', res)

			 allocate(this%node_bc(this%nnode_bc),stat=res)
			 if (res.ne.0)  rc = IOerr(lun, err_alloc, 'build_nodebc', &
					' temporary array temp', res)
			 m=1
			 do i=1,this%nedge_bc
				 do j=1,2
					 temp(m)=this%edge_iside(j,i)
					 m=m+1
				 end do
			 end do
			 call isort(2*this%nedge_bc,temp)
			 m=1
			 do i=1,2*this%nedge_bc,2
				 this%node_bc(m)=temp(i)
				 m=m+1
			 end do
			 deallocate(temp,stat=res)
			 if (res.ne.0)  rc = IOerr(lun, err_dealloc, 'build_nodebc', &
					' temporary array temp', res)
		 else
			 rc = IOerr(lun, err_val, 'build_nodebc', &
					' array edge_iside not built', res)
	    end if
	 else if (this%logical_dimension.eq.3) then
		 if (allocated(this%face_iside)) then
			 allocate(temp1(3*this%nface_bc),temp2(3*this%nface_bc),stat=res)
			 if (res.ne.0) rc = IOerr(lun, err_alloc, 'build_nodebc', &
				 'temporary array temp1 and temp2', res)
			 m = 1
			 do i = 1,this%nface_bc
				 do j = 1,3
					 temp1(m) = this%face_iside(j,i)
					 m = m + 1
				 end do
			 end do
			 call isort(3*this%nface_bc,temp1)
			 m = 1
			 temp2(1) = temp1(1)
			 do i = 2,3*this%nface_bc
				 if (temp1(i).ne.temp1(i-1)) then
					 m = m+1
					 temp2(m) = temp1(i)
				 end if
			 end do
			 this%nnode_bc = m
			 allocate(this%node_bc(this%nnode_bc),stat=res)
			 if (res.ne.0) rc = IOerr(lun, err_alloc, 'buil_nodebc', &
				 'array this%node_bc', res)
			 do i = 1,this%nnode_bc
				 this%node_bc(i) = temp2(i)
			 end do
			 deallocate(temp1,temp2,stat=res)
			 if (res.ne.0) rc = IOerr(lun, err_dealloc, 'build_nodebc', &
				 'temporary array temp1 and temp2', res)
		 else
			 rc = IOerr(lun, err_val, 'build_nodebc', &
					' array face_iside not built', res)
		 end if
    end if
    

  end subroutine build_nodebc

  !>-------------------------------------------------------------
  !> Build list of boundary nodes
  !> (private procedure for type mesh )
  !> 
  !> usage:
  !>     call 'var'%build_boundary_cell(lun)
  !>
  !> where:
  !> \param[in], optional, lun_err -> integer. I/O unit for error message output
  !<-------------------------------------------------------------
  subroutine build_boundary_cell(this,lun_err)
    use KindDeclaration
    use Globals
    implicit none
    class(abs_simplex_mesh), intent(inout) :: this
    integer, optional,       intent(in   ) :: lun_err
    !local
    logical :: rc
    integer :: res,lun
    integer :: i,j,m
    integer, allocatable :: temp(:)


    if (present(lun_err)) then
       lun=lun_err
    else
       lun=0
    end if

    if ( allocated(this%edge_plist)) then
       this%nnode_bc=this%nedge_bc
       allocate(temp(this%ncell),stat=res)
       if (res.ne.0)  rc = IOerr(lun, err_alloc, 'build_boundary_cell', &
            ' temporary array temp', res)

      
       this%nboundary_cell = 0
       temp=0
       do i=1,this%nedge_bc
          do j=1,2
             if ( this%edge_plist(j,i) .eq. 0 ) then
                temp(this%edge_plist(mod(j,2)+1,i))=1
             end if
          end do
       end do

       this%nboundary_cell=sum(temp)
       allocate(this%boundary_cell(this%nboundary_cell),stat=res)
       if (res.ne.0)  rc = IOerr(lun, err_alloc, 'build_boundary_cell', &
            ' temporary array temp', res)

       write(*,*)  this%nboundary_cell,this%ncell,size(this%boundary_cell)
       
       m=0
       do i=1,this%ncell
          if (temp(i) .ne. 0) then 
             m=m+1
             this%boundary_cell(m)=i
          end if
       end do

       
       deallocate(temp,stat=res)
       if (res.ne.0)  rc = IOerr(lun, err_dealloc, 'build_boundary_cell', &
            ' temporary array temp', res)
    else
       rc = IOerr(lun, err_val, 'build_boundary_cell', &
            ' array edge_iside not built', res)
    end if
    

  end subroutine build_boundary_cell


  
  !>-------------------------------------------------------------
  !> Procedure to reduce band width of matrix
  !> of node-node connection via Reverse Cuthill-Mckee ordering 
  !> (procedure private for type mesh)
  !> Build the topol and coord array of a conformally refined
  !> grid from a given coarser grid.
  !> Builds also the parent relationships.
  !>
  !> usage:
  !>     call 'var'%reorder_node([lun_err])
  !>
  !> where:
  !> \param[in ] (optional) lun_err -> integer. IO unit for err. msg. 
  !<-------------------------------------------------------------
  subroutine RCM_node_reorder(this,lun_err)
    use KindDeclaration
    use Globals
    use SparseMatrix
    implicit none
    class(abs_simplex_mesh), intent(inout) :: this
    integer, optional,       intent(in   ) :: lun_err
    ! local
    logical :: rc
    integer :: res,lun
    type(spmat) :: connection_matrix
    character(len=256) :: cell_type
    integer, allocatable :: perm(:),iperm(:)

    if (present(lun_err)) then
       lun=lun_err
    else
       lun=0
    end if
    !
    ! build node-node conection matrix
    !
    call this%build_nodenode_matrix(lun, .False.,connection_matrix)
    !
    ! generate ndoe permuatation
    !
    allocate(perm(this%nnode),iperm(this%nnode),stat=res)
    if (res.ne.0) rc = IOerr(lun, err_inp, 'data2grids', &
         ' work arrays perm iperm', res)
    call connection_matrix%genrcm(6,perm,iperm)
    call connection_matrix%kill(lun)

    !
    ! reorder nodes
    !
    call this%renumber(lun, this%nnode,perm,iperm)

    deallocate(perm,iperm,stat=res)
    if (res.ne.0) rc = IOerr(lun, err_dealloc, 'data2grids', &
         ' work arrays perm iperm', res)

  end subroutine RCM_node_reorder



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
  !> \param[in ] input_mesh -> type(mesh). mesh to be refined
  !>
  !<-------------------------------------------------------------
  subroutine build_refined_topol_coord(this, lun_err, input_mesh)
    use Globals
    implicit none
    ! External variables
    class(abs_simplex_mesh), intent(inout) :: this
    integer,     intent(in)  :: lun_err
    class(abs_simplex_mesh),  intent(in)  :: input_mesh
    ! Local vars
    type(abs_simplex_mesh) :: local_mesh 
    integer :: vert(4),edg(3),midnod(2),nod(3,2)
    integer :: iedge,inode,icell,ivert,izone,itetloc,ntri1
    integer :: res
    integer :: iedgeloc,jedgeloc,inodloc,i,j,jnode,ntetref(8),iloc,jloc
    logical :: rc
    integer :: n1,n2,n3,n4,n12,n13,n14,n23,n24,n34
    real(kind=double) :: dnrm2
    real(kind=double) :: axis_12_34(3),axis_13_24(3), axis_14_23(3),norms(3)
    integer :: imode
    select type ( input_mesh)
    type is (abs_simplex_mesh)
       local_mesh=input_mesh
    end select

    !
    ! free memory 
    !
    if ( this%is_initialized ) call this%kill(lun_err)

    

    if ( .not. allocated(local_mesh%edge_cnc) ) &
         call local_mesh%build_edge_connection(lun_err)
    
    ! assign subgrid dimension
    this%mesh_type    = local_mesh%mesh_type
    this%cell_type    = local_mesh%cell_type
    this%nnodeincell  = local_mesh%nnodeincell
    this%logical_dimension = local_mesh%logical_dimension
    this%ambient_dimension = local_mesh%ambient_dimension
    this%cell_id = local_mesh%cell_id


    this%nnode        = local_mesh%nnode + local_mesh%nedge
    this%ncell        = 2 ** (local_mesh%nnodeincell-1) * &
         local_mesh%ncell
    this%grid_level   = local_mesh%grid_level + 1
    this%nnode_parent = local_mesh%nnode
    this%ncell_parent = local_mesh%ncell

    ! allocate additional array for subgrid
    allocate (&
         this%topol(local_mesh%nnodeincell+1,this%ncell),&
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
       n1 = local_mesh%edge_iside(1,iedge)
       n2 = local_mesh%edge_iside(2,iedge)
       this%coord(:,local_mesh%nnode+iedge)=onehalf*(local_mesh%coord(:,n1) + local_mesh%coord(:,n2))
    end do

    ! subnode parent
    do inode = 1, local_mesh%nnode
       this%node_parent(1,inode) = inode  
       this%node_parent(2,inode) = inode
    end do
    do iedge = 1, local_mesh%nedge 
       inode = local_mesh%nnode + iedge
       this%node_parent(1,inode) = local_mesh%edge_iside(1,iedge)
       this%node_parent(2,inode) = local_mesh%edge_iside(2,iedge)
    end do


    select case  ( this%nnodeincell)
    case (3 ) 
       ! build new triang
       ntri1 = 0
       do icell=1,local_mesh%ncell
          ! store vertices
          do ivert = 1,3
             vert(ivert) = local_mesh%topol(ivert,icell)
          end do
          izone = local_mesh%topol(4,icell)
          ! store edges in triangle
          do iedgeloc=1,3
             edg(iedgeloc) = local_mesh%edge_cnc(iedgeloc,icell)
             ! store nodes in edges
             do inodloc = 1,2
                nod(iedgeloc,inodloc) = local_mesh%edge_iside(inodloc,edg(iedgeloc))
             end do
          end do
          ! build triangles
          do ivert = 1,3
             inode = vert(ivert)
             !find edges including inode
             iedgeloc = 0
             do jedgeloc = 1,3
                if ((inode .eq. nod(jedgeloc,1)) .or. &
                     (inode .eq. nod(jedgeloc,2))) then 
                   ! edg(jedgeloc) includes inod
                   iedgeloc = iedgeloc + 1
                   midnod(iedgeloc) = edg(jedgeloc)+local_mesh%nnode
                endif
             end do
             ! jedgeloc
             ntri1 = ntri1 + 1
             this%topol(1,ntri1) = inode
             do i = 1,2
                this%topol(i+1,ntri1) = midnod(i) 
             end do
             this%topol(4,ntri1) = izone
             this%cell_parent(ntri1) = icell
          end do                 ! iedgeloc
          ! build inner triangle
          ntri1 = ntri1 + 1
          do inodloc = 1,3
             iedgeloc = inodloc
             this%topol(inodloc,ntri1) = edg(iedgeloc)+local_mesh%nnode
          end do
          this%topol(4,ntri1) = izone
          this%cell_parent(ntri1) = icell
       end do

    case (4)
       ! initialize  ntetref
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
          ! build first four tetrahedrons
          this%topol(1,ntetref(1)) = n1 
          this%topol(2,ntetref(1)) = n12
          this%topol(3,ntetref(1)) = n13  
          this%topol(4,ntetref(1)) = n14 
          
          this%topol(1,ntetref(2)) = n2  
          this%topol(2,ntetref(2)) = n12
          this%topol(3,ntetref(2)) = n23
          this%topol(4,ntetref(2)) = n24 
          
          this%topol(1,ntetref(3)) = n3  
          this%topol(2,ntetref(3)) = n13 
          this%topol(3,ntetref(3)) = n23
          this%topol(4,ntetref(3)) = n34
          
          this%topol(1,ntetref(4)) = n4  
          this%topol(2,ntetref(4)) = n14 
          this%topol(3,ntetref(4)) = n24 
          this%topol(4,ntetref(4)) = n34

          ! evaluate Euclidean norm between opposite middle edge nodes
          ! and select the case with minimum distance
          axis_12_34=this%coord(:,n12)-this%coord(:,n34)
          axis_13_24=this%coord(:,n12)-this%coord(:,n34)
          axis_14_23=this%coord(:,n14)-this%coord(:,n23)
          norms(1)=dnrm2(3,axis_12_34,1)
          norms(2)=dnrm2(3,axis_13_24,1)
          norms(3)=dnrm2(3,axis_14_23,1)
          imode=1
          if (norms(2) < norms(1) ) imode=2
          if (norms(2) < norms(3) ) imode=3
          if (norms(3) < norms(1) ) imode=3
          if (norms(3) < norms(2) ) imode=3

          ! build other four tetrahedrons, whose construction
          ! depends on the norm evaluated above
          select case (imode)
          case (1)
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
          case(2)
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
          case(3)
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
          end select

          do iloc=1,8
             this%topol(5,ntetref(iloc)) = izone
             this%cell_parent(ntetref(iloc)) = icell
             ! sort the nodes of refined mesh
             ! in increasing order
             call isort(4,this%topol(1:4,ntetref(iloc)))
          end do
          

          ! update column numbers of the refined tetrahedrons
          ! in the array this%topol 
          ntetref(:) = ntetref(:) + 8

       end do
    end select
    
    this%is_initialized=.True.
    call local_mesh%kill(lun_err)

  end subroutine build_refined_topol_coord

  !>----------------------------------------------------------
  !> Procedure that reorders mesh nodes
  !> given a permutation vector.
  !> It works on members triang and coord
  !> and rebuilds all other member if necessary. 
  !>
  !> usage:
  !>     call 'var'%reorder(lun_err,perm_nodes,inv_perm_nodes)
  !>
  !> where:
  !> \param[in] lun_err   -> integer. I/O unit for error message output
  !> \param[in] n2permute -> integer. Length of permutation.
  !>                              n2permute=grid%nnode -> permute nodes
  !>                              n2permute=grid%cell -> permute cells
  !> \param[in] perm      -> integer. dimension(nnode)
  !>                                       node permutation array
  !> \param[in] inv_perm  -> integer. dimension(nnode)
  !>                                       inverse of perm_nodes
  !>------------------------------------------------------------------
  subroutine renumber_mesh(this,lun_err,n2permute,perm,inv_perm)
    use Globals
    implicit none
    class(abs_simplex_mesh), intent(inout) :: this
    integer,     intent(in   ) :: lun_err
    integer,     intent(in   ) :: n2permute
    integer,     intent(in   ) :: perm(n2permute) 
    integer,     intent(in   ) :: inv_perm(n2permute) 
    ! local 
    logical :: rc
    integer :: res
    integer :: inode,icell,i
    integer, allocatable :: topol_temp(:,:)
    integer ,allocatable :: new_node_parent(:,:)
    integer ,allocatable :: new_cell_parent(:)



    if ( n2permute .eq. this%nnode) then
       !  Permute the nodes according to the permutation vector.
       call double_col_permute(3, this%nnode, perm, this%coord)
       !
       !  Permute the node indices in the triangle array.
       do icell = 1, this%ncell
          do i = 1, this%nnodeincell
             inode = this%topol(i,icell)
             this%topol(i,icell) = inv_perm( inode )
          end do
       end do
       if ( this%grid_level> 0) then
          allocate(new_node_parent(2,this%nnode),stat=res)
          if (res .ne. 0) rc = IOerr(lun_err, err_alloc, 'renumber_mesh', &
               'work array new_node_parent')
          do inode=1,this%nnode
             do i=1,2
                new_node_parent(i,inv_perm(inode))=&
                     this%node_parent(i,inode)
             end do
          end do
          this%node_parent=new_node_parent
          deallocate(new_node_parent,stat=res)
          if (res .ne. 0) rc = IOerr(lun_err, err_dealloc, 'renumber_mesh', &
               'work array new_node_parent')
       end if

    else  if ( n2permute .eq. this%ncell) then
       allocate(topol_temp(this%nnodeincell+1,this%ncell),stat=res)
       if (res.ne.0) rc = IOerr(lun_err, err_val, 'renumber_mesh', &
            ' temp array topol_temp')
       topol_temp=this%topol
       do icell=1,this%ncell
          this%topol(:,icell) = topol_temp(:,perm(icell))
       end do
       deallocate(topol_temp,stat=res)
       if (res.ne.0) rc = IOerr(lun_err, err_dealloc, 'renumber_mesh', &
            ' temp array topol_temp')

       if ( this%grid_level> 0) then
          allocate(new_cell_parent(this%ncell)) 
          do icell=1,this%nnode
             new_cell_parent(icell)=perm(this%cell_parent(icell))
          end do
          this%cell_parent=new_cell_parent
          deallocate(new_cell_parent) 
       end if


    else
       rc = IOerr(lun_err, err_val, 'renumber_mesh', &
            'permutation length not equal to grid nodes or cell number')
    end if





  end subroutine renumber_mesh


  !>-----------------------------------------------------
  !> Minimal mesh parameter
  !>-----------------------------------------------------
  function meshpar(this,flag) result(h)
    use Globals
    implicit none
    class(abs_simplex_mesh), intent(in) :: this
    integer,     intent(in) :: flag 
    real(kind=double) :: h
    !local

    select case (this%nnodeincell)
    case (3)
       if (flag .eq. 0) h=sqrt(2*sum(this%size_cell)/this%ncell)
       if (flag .eq. 1) h=minval(sqrt(this%size_cell)/2.0d0)
       if (flag .eq. 2) h=maxval(sqrt(this%size_cell)/2.0d0)
       if (flag .eq. 3) h=sum(this%size_edge)*one/this%nedge
       if (flag .eq. 4) h=minval(this%size_edge)
       if (flag .eq. 5) h=maxval(this%size_edge)
    case default
       write(*,*) 'MeshPar defined only for tringulation'
    end select
  end function meshpar


  !>-------------------------------------------------------------
  !> Part of the static constructor.
  !> (procedure public for type abs_simplex_mesh)
  !> Instantiate and initialize type abs_simplex_mesh
  !> from another via a selector of nodes or cells
  !>
  !> usage:
  !>     call this%read_mesh(lun_err, file2read)
  !>
  !> where:
  !> \param[in] lun_err     -> integer. IO unit for error msg
  !> \param[in] parent_mesh -> type(abs_mesh). Original mesh
  !> \param[in] ndata       -> integer. Lenght of selector
  !>                           Must be equal to ncel or nnode
  !> \param[in] selector    -> logical(dim =ndata) 
  !>                           If selector(i)=.False. remove node(i)/cell(i) 
  !<-------------------------------------------------------------
  subroutine init_from_data(this, lun_err, &
       nnode,ncell, nnodeincell,cell_type,&
       topol, coord,logical_dimension)
    use Globals
    class(abs_simplex_mesh), intent(inout) :: this
    integer,                  intent(in   ) :: lun_err
    integer,                  intent(in   ) :: nnode
    integer,                  intent(in   ) :: ncell
    integer,                  intent(in   ) :: nnodeincell
    character(len=*),         intent(in   ) :: cell_type
    integer,                  intent(in   ) :: topol(nnodeincell,ncell)
    real(kind=double),        intent(in   ) :: coord(3,nnode)
    integer,optional,         intent(in   ) :: logical_dimension
    

    ! local
    logical :: rc
    integer :: res
    integer :: icell
    real(kind=double) :: zmin, zmax
    character(len=256) :: clean

    if ( this%is_initialized ) call this%kill(lun_err)
    
    this%nnode       = nnode
    this%ncell       = ncell
    this%nnodeincell = nnodeincell
    this%cell_type   = etb(cell_type)

    clean=etb(cell_type)

    zmin=minval(coord(3,:))
    zmax=maxval(coord(3,:))
    if (abs(zmin-zmax)<1e-15 ) then
       this%ambient_dimension=2           
    else
       this%ambient_dimension=3
    end if


    select case (clean) 
    case( 'triangle')      
       if (this%ambient_dimension .eq. 2  ) then
          this%logical_dimension=2
          this%mesh_type='2d'
       else
          this%logical_dimension=3
          if ( present(logical_dimension)) this%logical_dimension=logical_dimension
          this%mesh_type = 'surface'
       end if
       this%cell_id = 5
    case ('tetrahedron')
       this%logical_dimension=3
       this%cell_id = 10
       this%mesh_type='3d'
    case ('square')
       this%logical_dimension=2
       this%cell_id = 9
       this%mesh_type='squares'
    case ('edge')
       this%logical_dimension=1
       this%cell_id = 3
       this%mesh_type='graph'
    end select

    allocate(this%topol(nnodeincell+1,ncell),&
         this%coord(3,nnode),stat=res)
    if(res .ne. 0) rc = IOerr(lun_err, err_alloc, 'init_from_data', &
         'type abs_simplex_mesh member topol coord')

    do icell=1,ncell
       this%topol(1:nnodeincell,icell)=topol(1:nnodeincell,icell)
       this%topol(1+nnodeincell,icell)=icell
    end do

    this%coord = coord

    this%is_initialized=.True.

  end subroutine init_from_data

  subroutine mesh_constructor_from_data(this, lun_err, &
       nnode,ncell, nnodeincell,cell_type,&
       topol, coord)

    type(abs_simplex_mesh),   intent(inout) :: this
    integer,                  intent(in   ) :: lun_err
    integer,                  intent(in   ) :: nnode
    integer,                  intent(in   ) :: ncell
    integer,                  intent(in   ) :: nnodeincell
    character(len=*),         intent(in   ) :: cell_type
    integer,                  intent(in   ) :: topol(nnodeincell,ncell)
    real(kind=double),        intent(in   ) :: coord(3,nnode)

    call this%init_from_data(lun_err, &
         nnode,ncell, nnodeincell,cell_type,&
         topol, coord)
  end subroutine mesh_constructor_from_data


  subroutine eval_node2cell(this,nodevar,cellvar)
    use Globals
    implicit none
    ! vars
    class(abs_simplex_mesh), intent(in) :: this
    real(kind=double), intent(in) :: nodevar(this%nnode)
    real(kind=double), intent(out) :: cellvar(this%ncell)
    ! local vars
    integer :: icell, inode, iloc

    cellvar = zero
    do icell=1,this%ncell

       do iloc=1,this%nnodeincell
          inode = this%topol(iloc,icell)

          cellvar(icell) = cellvar(icell) &
               + nodevar(inode) / this%nnodeincell
       end do
    end do

  end subroutine eval_node2cell

  subroutine eval_cell2node(this,cellvar,nodevar)
    use Globals
    implicit none
    ! vars
    class(abs_simplex_mesh), intent(in) :: this
    real(kind=double), intent(in) :: cellvar(this%ncell)
    real(kind=double), intent(out) :: nodevar(this%nnode)
    ! local vars
    integer :: icell, inode, iloc
    real(kind=double) :: loc_val

    nodevar = zero
    do icell=1,this%ncell
       loc_val = cellvar(icell)*this%size_cell(icell)/this%nnodeincell

       do iloc=1,this%nnodeincell
          inode = this%topol(iloc,icell)

          nodevar(inode) = nodevar(inode) + loc_val
       end do
    end do
    do inode=1,this%nnode
       nodevar(inode) = nodevar(inode)/this%size_node(inode)
    end do
    
  end subroutine eval_cell2node

  subroutine eval_zone2node(this,zonevar,nodevar)
    use Globals
    implicit none
    ! vars
    class(abs_simplex_mesh), intent(in) :: this
    real(kind=double), intent(in) :: zonevar(this%nzone)
    real(kind=double), intent(inout) :: nodevar(this%nnode)
    ! local vars
    integer :: icell, inode, iloc
    real(kind=double) :: loc_val

    nodevar = zero
    do icell=1,this%ncell
       loc_val = zonevar(this%topol(this%nnodeincell+1,icell))
       loc_val = loc_val*this%size_cell(icell)/this%nnodeincell

       do iloc=1,this%nnodeincell
          inode = this%topol(iloc,icell)

          nodevar(inode) = nodevar(inode) + loc_val
       end do
    end do
    do inode=1,this%nnode
       nodevar(inode) = nodevar(inode)/this%size_node(inode)
    end do
    
  end subroutine eval_zone2node




    subroutine circumcenter(v1,v2,v3,c)
    use Globals
    implicit none
    real(kind=double),  intent(in ) :: v1(2),v2(2),v3(2)
    real(kind=double),  intent(out) :: c(2)
    !local
    real(kind=double) :: d
    
    d=2.d0*(&
          v1(1)*( v2(2)-v3(2) ) + &
          v2(1)*( v3(2)-v1(2) ) + &
          v3(1)*( v1(2)-v2(2) ) )
 
    c(1)=(&
         (v1(1)**2+v1(2)**2) * (v2(2)-v3(2)) + &
         (v2(1)**2+v2(2)**2) * (v3(2)-v1(2)) + &
         (v3(1)**2+v3(2)**2) * (v1(2)-v2(2))  )  /d
    
    c(2)=(&
         (v1(1)**2+v1(2)**2) * (v3(1)-v2(1)) + &
         (v2(1)**2+v2(2)**2) * (v1(1)-v3(1)) + &
         (v3(1)**2+v3(2)**2) * (v2(1)-v1(1))  )  /d
    
  end subroutine circumcenter

  function scalar_product(this,vector1,vector2) result(out)
    use Globals
    implicit none
    class(abs_simplex_mesh), intent(in   ) :: this
    real(kind=double),       intent(in   ) :: vector1(this%logical_dimension)
    real(kind=double),       intent(in   ) :: vector2(this%logical_dimension)
    real(kind=double) :: out
    !
    real(kind=double) :: ddot
    
    out=ddot(this%logical_dimension,vector1,1,vector2,1)

  end function scalar_product

  ! function getting the number of nodes for a given cell
  function get_nnodes_in_cell(this, icell) result(nnodeincell)
    implicit none
    class(abs_simplex_mesh), intent(in   ) :: this
    integer,                 intent(in   ) :: icell
    integer :: nnodeincell    
    nnodeincell=this%nnodeincell

  end function get_nnodes_in_cell

  ! get the list of nodes in a given cell
  subroutine get_nodes_in_cell(this, icell, nnodeincell, nodes)
    implicit none
    class(abs_simplex_mesh), intent(in   ) :: this
    integer,                 intent(in   ) :: icell
    integer,                 intent(in   ) :: nnodeincell
    integer,                 intent(inout) :: nodes(nnodeincell)
    !local
    integer :: i
     
    do i=1,nnodeincell
       nodes(i) = this%topol(i,icell)
    end do

  end subroutine get_nodes_in_cell

  ! get the list of nodes in a given cell
  subroutine get_coord_of_node(this, inode, coord)
    implicit none
    class(abs_simplex_mesh), intent(in   ) :: this
    integer,                 intent(in   ) :: inode
    real(kind=double),       intent(inout) :: coord(3)
    ! local
    integer :: i

    do i = 1,this%ambient_dimension
       coord(i) = this%coord(i,inode)
    end do

  end subroutine get_coord_of_node

    

  function metric_scalar_product(this,vector1,vector2) result(out)
    use Globals
    implicit none
    class(metric_mesh), intent(in   ) :: this
    real(kind=double),       intent(in   ) :: vector1(this%logical_dimension)
    real(kind=double),       intent(in   ) :: vector2(this%logical_dimension)
    real(kind=double) :: out
    !
    real(kind=double) :: ddot
    real(kind=double) :: temp(2)

    temp=matmul(this%metric,vector2)
    
    out=ddot(this%logical_dimension,vector1,1,temp,1)

  end function metric_scalar_product

  
!!$  subroutine build_cartesian2barycentric(this,lun_err,info,&
!!$       icell,inverse_trasform,trasform)
!!$    use Globals
!!$    implicit none
!!$    class(abs_simplex_mesh), intent(in   ) :: this
!!$    integer,                 intent(in   ) :: lun_err
!!$    integer,                 intent(inout) :: info
!!$    integer,                 intent(in   ) :: icell
!!$    type(densemat),          intent(inout) :: inverse_trasform
!!$    type(denseprec),         intent(inout) :: trasform
!!$    !
!!$    real(kind=double) :: ddot
!!$    real(kind=double) :: temp(2)
!!$    integer :: i,j
!!$    
!!$    ! init tranform from 
!!$    call inverse_trasform%init(lun_err,this%nnodeincell,this%nnodeincell)
!!$
!!$    this%nnodeincell
!!$    do i=1, this%nnodeincell
!!$       inverse_trasform%coeff(1,i)=one
!!$       do j=2, this%nnodeincell
!!$          inode=
!!$          inverse_trasform%coeff(j,i)=this%coord(1,inode)
!!$       end do
!!$       
!!$    end  do
!!$  end subroutine build_cartesian2barycentric

  !> Solution taken from
  !> https://answers.unity.com/questions/383804/calculate-uv-coordinates-of-3d-point-on-plane-of-m.html

  subroutine cartesian2barycentric(this,icell,cartesian, barycentric)
    use Globals
    implicit none
    class(abs_simplex_mesh), intent(in   ) :: this
    integer,                 intent(in   ) :: icell
    real(kind=double),       intent(inout) :: cartesian(3)
    real(kind=double),       intent(inout) :: barycentric(this%nnodeincell)
    !
    integer :: i
    real(kind=double) :: nodes_coordinate(3,4) ! 4 is to handle 3d data

    do i=1,this%nnodeincell
       nodes_coordinate(:,i)=this%coord(:,this%topol(i,icell))
    end do
    call pure_cartesian2barycentric(&
       this%nnodeincell,&
       nodes_coordinate,&
       cartesian,&
       barycentric)
    
  end subroutine cartesian2barycentric

  subroutine pure_cartesian2barycentric(&
       nnodeincell,&
       nodes_coordinate,&
       cartesian,&
       barycentric)
    implicit none
    integer :: nnodeincell
    real(kind=double) :: nodes_coordinate(3,4)
    real(kind=double) :: cartesian(3)
    real(kind=double) :: barycentric(nnodeincell)
    ! local 
    real(kind=double) :: f(3,3),v(3,3),temp1(3),temp2(3),va(3),area
    real(kind=double) :: dnrm2,ddot
    integer :: i,j

    ! compute vectors v_i - p
    do i=1,nnodeincell
       f(:,i)=nodes_coordinate(:,i)-cartesian(:)
    end do
    ! compute vector ortogonal to the triangle
    temp1 = nodes_coordinate(:,1)-nodes_coordinate(:,2)
    temp2 = nodes_coordinate(:,1)-nodes_coordinate(:,3)
    va    = cross(temp1,temp2)
    area  = dnrm2(3,va,1)

    ! compute cross 
    v(:,1)=cross(f(:,2),f(:,3))
    v(:,2)=cross(f(:,3),f(:,1))
    v(:,3)=cross(f(:,1),f(:,2))

    ! compute singed weights
    do i=1,nnodeincell
       barycentric(i)=dnrm2(3,v(:,i),1)/area*&
            sign(one,ddot(3,va,1,v(:,i),1))
    end do
  
  end subroutine pure_cartesian2barycentric
  
  
end module AbstractGeometry
