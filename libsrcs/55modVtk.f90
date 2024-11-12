!>-------------------------------------------------------------
!> procedures for writing vtk files. 
!>
!> only mesh and scalars up to now
!<------------------------------------------------------------- 
module vtkloc
  use Globals
  implicit none
  
  private
  !> outputs mesh info in vtk format
  public :: vtkmesh
  !> outputs scalar vars in vtk format
  public :: vtkscal
  !> outputs scalar vars in vtk format
  public :: vtkvectors
  !> inputs scalar vars from vtk format
  public :: scalin
!!$  !> outputs mesh info in vtk format
!!$  public :: vtkgraph
contains

  !>-------------------------------------------------------------
  !> outputs mesh info in vtk format
  !>
  !> usage:
  !>     call vtkmesh(lun,fname,fileformat,time,geo)
  !>
  !> where:
  !> \param[in] lun: (integer) unit for error message output
  !> \param[in] fname: (character) filename for main output
  !> \param[in] fileformat: (character) 'ASCII' or 'BINARY'
  !> \param[in] time: (double) simulation time of output
  !> \param[in] geo: (mesh) mesh structure
  !>
  !> uses private subroutines alloc_mesh_tmp and
  !>      dealloc_mesh_tmp
  !<-----------------------------------------------------------
  subroutine vtkmesh(lun,fname,fileformat,time,geo)
    use Globals
    use AbstractGeometry
    use LIB_VTK_IO
    implicit none
    ! vars
    integer, intent(in) :: lun
    character(len=*),  intent(in) :: fname,fileformat
    real(kind=single), intent(in) :: time
    type(abs_simplex_mesh),        intent(in) :: geo
    ! local vars
    real(kind=double), allocatable :: xcoord(:),ycoord(:),zcoord(:)
    integer, allocatable :: cells(:),ctype(:),material(:)
    integer :: ncell,nnod,nodes_in_tria
    integer :: res,E_IO
    integer :: i
    logical :: rc
    
    ncell = geo%ncell
    nnod = geo%nnode
    nodes_in_tria = geo%nnodeincell
    call alloc_mesh_tmp(lun,ncell,nnod,nodes_in_tria,geo,&
         cells,ctype,xcoord,ycoord,zcoord)

    allocate(material(ncell),stat=res)
    if(res .ne. 0) rc = IOerr(lun, err_alloc, 'vtkmesh', &
         '  material (local array)',res)

    do i=1,ncell
       material(i) = geo%topol(nodes_in_tria+1,i)
    end do

    E_IO = VTK_INI(output_format = fileformat, &
         filename = trim(fname), &
         title = 'Unstructured Grid', &
         time = time, &
         mesh_topology = 'UNSTRUCTURED_GRID')
    if (E_IO .ne. 0) rc = IOerr(lun, err_vtk, 'vtkmesh', &
         ' error in initializing vtk file, E_IO=',E_IO)
    E_IO = VTK_GEO(Nn = nnod, X = xcoord, Y = ycoord, Z = zcoord)
    if (E_IO .ne. 0) rc = IOerr(lun, err_vtk, 'vtkmesh', &
         ' error in writing coordinates in vtk file, E_IO=',E_IO)
    E_IO = VTK_CON(NC = ncell, connect = cells, cell_type = ctype)
    if (E_IO .ne. 0) rc = IOerr(lun, err_vtk, 'vtkmesh', &
         ' error in writing mesh topology in vtk file, E_IO=',E_IO)
!!$    ! Materiali
!!$    E_IO = VTK_DAT(NC_NN = ncell, var_location = 'cell')
!!$    if (E_IO .ne. 0) rc = IOerr(lun, err_vtk, 'vtkmesh', &
!!$         ' error in writing material headers in vtk file, E_IO=',E_IO)
!!$    E_IO = VTK_VAR(NC_NN = ncell, varname = 'Material', var = material)
!!$    if (E_IO .ne. 0) rc = IOerr(lun, err_vtk, 'vtkmesh', &
!!$         ' error in writing mesh materials in vtk file, E_IO=',E_IO)

    E_IO = VTK_END()
    if (E_IO .ne. 0) rc = IOerr(lun, err_vtk, 'vtkmesh', &
         ' error in finalizing vtk file, E_IO=',E_IO)

    call dealloc_mesh_tmp(lun,cells,ctype,xcoord,ycoord,zcoord) 

    deallocate(material,stat=res)
    if (res.ne.0) rc=IOerr(lun, err_dealloc, 'vtkmesh', &
         'material (local array)',res)
  end subroutine vtkmesh

  !>-------------------------------------------------------------
  !> outputs scalar var in vtk format
  !>
  !> usage:
  !>     call vtkscal(lun,fname,fileformat,time,geo,
  !>                  nvar,utype,uname,uvar)
  !>
  !> where:
  !> \param[in] lun: (integer) unit for error message output
  !> \param[in] fname: (character) filename for main output
  !> \param[in] fileformat: (character) 'ASCII' or 'BINARY'
  !> \param[in] time: (double) simulation time of output
  !> \param[in] geo: (mesh) mesh structure
  !> \param[in] nvar: (integer) number of vars
  !> \param[in] mazsize: (integer) max sazie of the variables
  !> \param[in] utype: (character array, dim nvar) var types
  !> \param[in] uname: (character array, dim nvar) var names
  !> \param[in] uvar: (double array, dim (nvar,:)) variable
  !>
  !> uses private subroutines alloc_mesh_tmp and
  !>      dealloc_mesh_tmp
  !<-----------------------------------------------------------
  subroutine vtkscal(lun,fname,fileformat,time,geo, &
       nvar,ntot,utype,uname,uvar)
    use Globals
    use AbstractGeometry
    use LIB_VTK_IO
    implicit none
    ! vars
    integer :: lun
    character(len=*), intent(in) :: fname,fileformat
    real(kind=single), intent(in) :: time
    type(abs_simplex_mesh), intent(in) :: geo
    integer, intent(in) :: nvar,ntot
    character(len=*), intent(in) :: utype(nvar),uname(nvar)
    real(kind=double), intent(in) :: uvar(ntot)
    ! local vars
    real(kind=double), allocatable :: xcoord(:),ycoord(:),zcoord(:),uuu(:)
    integer, allocatable :: cells(:),ctype(:)    
    character(len=256) :: msg
    character(len=10) :: rdwr
    character(len=4) :: case_arg
    integer :: ncell,nnod,nodes_in_tria
    integer :: res,E_IO
    integer :: i
    integer :: size_array(nvar),start(nvar),finish(nvar)
    logical :: var_cell,var_node,rc
    
    ncell = geo%ncell
    nnod = geo%nnode
    nodes_in_tria = geo%nnodeincell

    
    call alloc_mesh_tmp(lun,ncell,nnod,nodes_in_tria,geo,&
         cells,ctype,xcoord,ycoord,zcoord)
    E_IO = VTK_INI(output_format = fileformat, &
         filename = trim(fname), &
         title = 'Unstructured Grid', &
         time = time, &
         mesh_topology = 'UNSTRUCTURED_GRID')
    if (E_IO .ne. 0) rc = IOerr(lun, err_vtk, 'vtkscal', &
         ' error in initializing vtk file, E_IO=',E_IO)
    E_IO = VTK_GEO(Nn = nnod, X = xcoord, Y = ycoord, Z = zcoord)
    if (E_IO .ne. 0) rc = IOerr(lun, err_vtk, 'vtkscal', &
         ' error in writing coordinates in vtk file, E_IO=',E_IO)
    E_IO = VTK_CON(NC = ncell, connect = cells, cell_type = ctype)
    if (E_IO .ne. 0) rc = IOerr(lun, err_vtk, 'vtkscal', &
         ' error in writing mesh topology in vtk file, E_IO=',E_IO)
    
    var_cell=.false.
    var_node=.false.
    size_array=0
    start=1
    finish=1
    do i=1,nvar
       case_arg = trim(adjustl(utype(i)))
       if (case_arg(1:4).eq.'cell') then
          var_cell=.true.
          size_array(i)=ncell
       end if
       if (case_arg(1:4).eq.'node') then
          var_node=.true.
          size_array(i)=nnod
       end if

       if ( i .ne. 1) start(i) = finish(i-1) + 1
       finish(i)=start(i)+size_array(i)-1       
       !write(*,*) case_arg, size_array(i),start(i),finish(i) 
    end do
   

    if(var_cell) then
       ! write out cell data types (vectors or scalars)
       E_IO = VTK_DAT(NC_NN = ncell, var_location = 'cell')
       if (E_IO .ne. 0) then
          msg = 'writing headers for var ='
          write(rdwr,'(a)') uname(i)
          msg = trim(msg) // trim(rdwr) // ' at component '
          write(rdwr,'(i5)') i
          msg = trim(msg) // trim(rdwr) // ' of '
          write(rdwr,'(i5)') nvar
          msg = trim(msg) // trim(rdwr) // '. E_IO= '
          write(rdwr,'(i5)') E_IO
          msg = trim(msg) // trim(rdwr) 
          rc = IOerr(lun, err_vtk, 'vtkscal', msg, E_IO)
       end if
       allocate(uuu(ncell),stat=res)
       if(res .ne. 0) rc = IOerr(lun, err_alloc, 'vtkscal', &
            'uuu (local array)',res)
       do i=1,nvar
          case_arg = trim(adjustl(utype(i)))
          if(case_arg.eq.'cell') then
!!$             do j=1,ncell
!!$                uuu(j)=uvar(j,i)
!!$             end do
             E_IO = VTK_VAR(NC_NN = ncell, varname = trim(uname(i)),&
                  var = uvar(start(i):finish(i)) )
             if (E_IO .ne. 0) then
                msg = 'writing headers for var ='
                write(rdwr,'(a)') uname(i)
                msg = trim(msg) // trim(rdwr) // ' at component '
                write(rdwr,'(i5)') i
                msg = trim(msg) // trim(rdwr) // ' of '
                write(rdwr,'(i5)') nvar
                msg = trim(msg) // trim(rdwr) // '. E_IO= '
                write(rdwr,'(i5)') E_IO
                msg = trim(msg) // trim(rdwr) 
                rc = IOerr(lun, err_vtk, 'vtkscal', msg, E_IO)
             end if
          end if
       end do
       deallocate(uuu,stat=res)
       if (res.ne.0) rc=IOerr(lun, err_dealloc, 'vtkscal', &
            'uuu (local array)', res)
    end if
    if(var_node) then
       ! write out node data types (vectors or scalars)
       E_IO = VTK_DAT(NC_NN = nnod, var_location = 'node')
       if (E_IO .ne. 0) then
          msg = 'writing headers for var ='
          write(rdwr,'(a)') uname(i)
          msg = trim(msg) // trim(rdwr) // ' at component '
          write(rdwr,'(i5)') i
          msg = trim(msg) // trim(rdwr) // ' of '
          write(rdwr,'(i5)') nvar
          msg = trim(msg) // trim(rdwr) // '. E_IO= '
          write(rdwr,'(i5)') E_IO
          msg = trim(msg) // trim(rdwr) 
          rc = IOerr(lun, err_vtk, 'vtkscal', msg, E_IO)
       end if
       allocate(uuu(nnod),stat=res)
       if(res .ne. 0) rc = IOerr(lun, err_alloc, 'vtkscal', &
            'uuu (local array)',res)
       do i=1,nvar
          case_arg = trim(adjustl(utype(i)))
          if(case_arg.eq.'node') then
!!$             do j=1,nnod
!!$                uuu(j)=uvar(j,i)
!!$             end do
             E_IO = VTK_VAR(NC_NN = nnod, varname = trim(uname(i)),&
                  var = uvar(start(i):finish(i)) )
             if (E_IO .ne. 0) then
                msg = 'writing headers for var ='
                write(rdwr,'(a)') uname(i)
                msg = trim(msg) // trim(rdwr) // ' at component '
                write(rdwr,'(i5)') i
                msg = trim(msg) // trim(rdwr) // ' of '
                write(rdwr,'(i5)') nvar
                msg = trim(msg) // trim(rdwr) // '. E_IO= '
                write(rdwr,'(i5)') E_IO
                msg = trim(msg) // trim(rdwr) 
                rc = IOerr(lun, err_vtk, 'vtkscal', msg, E_IO)
             end if
          end if
       end do
       deallocate(uuu,stat=res)
       if (res.ne.0) rc=IOerr(lun, err_dealloc, 'vtkscal', &
            'uuu (local array)')
    end if
    
    E_IO = VTK_END()
    if (E_IO .ne. 0) then
       write(rdwr,'(i5)') E_IO
       msg = 'error in finalizing vtk file in vtkscal, E_IO='
       msg = trim(msg) // trim(rdwr)
       rc = IOerr(lun, err_vtk, 'vtkscal', msg, E_IO)
    end if
    
    call dealloc_mesh_tmp(lun,cells,ctype,xcoord,ycoord,zcoord) 
    
    return
  end subroutine vtkscal

  !>-------------------------------------------------------------
  !> outputs scalar var in vtk format
  !>
  !> usage:
  !>     call vtkscal(lun,fname,fileformat,time,geo,
  !>                  nvar,utype,uname,uvar)
  !>
  !> where:
  !> \param[in] lun: (integer) unit for error message output
  !> \param[in] fname: (character) filename for main output
  !> \param[in] fileformat: (character) 'ASCII' or 'BINARY'
  !> \param[in] time: (double) simulation time of output
  !> \param[in] geo: (mesh) mesh structure
  !> \param[in] nvar: (integer) number of vars
  !> \param[in] mazsize: (integer) max sazie of the variables
  !> \param[in] utype: (character array, dim nvar) var types
  !> \param[in] uname: (character array, dim nvar) var names
  !> \param[in] uvar: (double array, dim (nvar,:)) variable
  !>
  !> uses private subroutines alloc_mesh_tmp and
  !>      dealloc_mesh_tmp
  !<-----------------------------------------------------------
  subroutine vtkvectors(lun,fname,fileformat,time,geo, &
       ndata,utype,uname,xvar,yvar,zvar)
    use Globals
    use AbstractGeometry
    use LIB_VTK_IO
    implicit none
    ! vars
    integer :: lun
    character(len=*), intent(in) :: fname,fileformat
    real(kind=single), intent(in) :: time
    type(abs_simplex_mesh), intent(in) :: geo
    integer, intent(in) :: ndata
    character(len=*), intent(in) :: utype,uname
    real(kind=double), intent(in) :: xvar(ndata),yvar(ndata),zvar(ndata)
    ! local vars
    real(kind=double), allocatable :: xcoord(:),ycoord(:),zcoord(:),uuu(:)
    integer, allocatable :: cells(:),ctype(:)    
    character(len=256) :: msg
    character(len=256) :: uname_copy
    character(len=10) :: rdwr
    character(len=4) :: case_arg
    integer :: ncell,nnod,nodes_in_tria
    integer :: res,E_IO
    integer :: i
    logical :: var_cell,var_node,rc
    
    ncell = geo%ncell
    nnod = geo%nnode
    nodes_in_tria = geo%nnodeincell
    
    call alloc_mesh_tmp(lun,ncell,nnod,nodes_in_tria,geo,&
         cells,ctype,xcoord,ycoord,zcoord)
    E_IO = VTK_INI(output_format = fileformat, &
         filename = trim(fname), &
         title = 'Unstructured Grid', &
         time = time, &
         mesh_topology = 'UNSTRUCTURED_GRID')
    if (E_IO .ne. 0) rc = IOerr(lun, err_vtk, 'vtkscal', &
         ' error in initializing vtk file, E_IO=',E_IO)
    E_IO = VTK_GEO(Nn = nnod, X = xcoord, Y = ycoord, Z = zcoord)
    if (E_IO .ne. 0) rc = IOerr(lun, err_vtk, 'vtkscal', &
         ' error in writing coordinates in vtk file, E_IO=',E_IO)
    E_IO = VTK_CON(NC = ncell, connect = cells, cell_type = ctype)
    if (E_IO .ne. 0) rc = IOerr(lun, err_vtk, 'vtkscal', &
         ' error in writing mesh topology in vtk file, E_IO=',E_IO)
    
    ! write out node data types (vectors or scalars)
    !uname_copy=etb(uname)
    E_IO = VTK_DAT(NC_NN = ndata, var_location = utype)
!!$ varname = uname_copy, &
!!$         varX=xvar, varY=yvar, varZ=zvar)
!!$    E_IO = VTK_DAT(NC_NN = ndata, utype) varname = uname_copy, &
!!$         varX=xvar, varY=yvar, varZ=zvar)

    if (E_IO .ne. 0) then
       msg = 'writing headers for var ='
       write(rdwr,'(a)') uname
       msg = trim(msg) // trim(rdwr) // ' at component '
       write(rdwr,'(i5)') i
       msg = trim(msg) // trim(rdwr) // ' of '
       write(rdwr,'(i5)') 1
       msg = trim(msg) // trim(rdwr) // '. E_IO= '
       write(rdwr,'(i5)') E_IO
       msg = trim(msg) // trim(rdwr) 
       rc = IOerr(lun, err_vtk, 'vtkscal', msg, E_IO)
    end if
    E_IO = VTK_VAR('vect', NC_NN = ndata, varname = uname,&
        varX=xvar, varY=yvar, varZ=zvar)

       
    E_IO = VTK_END()
    if (E_IO .ne. 0) then
       write(rdwr,'(i5)') E_IO
       msg = 'error in finalizing vtk file in vtkscal, E_IO='
       msg = trim(msg) // trim(rdwr)
       rc = IOerr(lun, err_vtk, 'vtkscal', msg, E_IO)
    end if
    
    call dealloc_mesh_tmp(lun,cells,ctype,xcoord,ycoord,zcoord) 
    
    return
  end subroutine vtkvectors



  !>-------------------------------------------------------------
  !> outputs scalar var in vtk format
  !>
  !> usage:
  !>     call vtkscal(lun,fname,fileformat,time,geo,
  !>                  nvar,utype,uname,uvar)
  !>
  !> where:
  !> \param[in] lun: (integer) unit for error message output
  !> \param[in] fname: (character) filename for main output
  !> \param[in] fileformat: (character) 'ASCII' or 'BINARY'
  !> \param[in] time: (double) simulation time of output
  !> \param[in] geo: (mesh) mesh structure
  !> \param[in] nvar: (integer) number of vars
  !> \param[in] mazsize: (integer) max sazie of the variables
  !> \param[in] utype: (character array, dim nvar) var types
  !> \param[in] uname: (character array, dim nvar) var names
  !> \param[in] uvar: (double array, dim (nvar,:)) variable
  !>
  !> uses private subroutines alloc_mesh_tmp and
  !>      dealloc_mesh_tmp
  !<-----------------------------------------------------------
  subroutine scalin(lun_err,lun_vtk,&
       ndata,&
       fname,&
       name_data,&
       data)
    use Globals
    implicit none
    ! vars
    integer,           intent(in ) :: lun_err
    integer,           intent(in ) :: lun_vtk
    integer,           intent(in ) :: ndata
    character(len=*),  intent(in ) :: fname
    character(len=*),  intent(in ) :: name_data
    real(kind=double), intent(out) :: data(ndata)

    ! local vars
    logical :: rc,found_DATA
    character(len=256) :: input_text ,scratch
    character(len=256) :: text_DATA,clean_name
    character(len=7)   :: SCALARS
    !character, allocatable :: text_DATA(:), clean_name(:)
    integer :: i,res

    
    open(lun_vtk,file = etb(fname),iostat = res)
    if(res .ne. 0) rc = IOerr(lun_err, err_inp , 'vtk2data', &
         fname // ' local var text_DATA, clean_name',res)
    
    clean_name=etb(name_data)    

    found_DATA = .false.

    ! found the line where to start and check dimension
    do while (.not. found_DATA)
       read(lun_vtk,'(a256)',iostat=res) input_text
       if(res .ne. 0) rc = IOerr(lun_err, err_inp , 'scalin', &
         (fname // ' local var input_text'),res)
       !write(*,*) etb(input_text)
       !read first word
       !read(lun_vtk,'(a256)',iostat=res) input_text
       read(input_text,*) SCALARS
       if ( SCALARS .EQ. 'SCALARS' ) then
          read(input_text,*) SCALARS, text_DATA 
          if ( etb(text_DATA) .eq. clean_name ) then
             found_DATA = .true.
             write(*,*) 'found variable "', etb(text_DATA),'"'
             read(lun_vtk,'(a256)',iostat=res) scratch
          end if
       end if
    end do
    
    do i = 1,ndata
       read( lun_vtk,*,iostat=res) data(i)
       if(res .ne. 0) THEN
          rc = IOerr(lun_err, err_inp , 'scalin', &
               trim(fname) // &
               ' real array data',res)
       end if
    end do
    close (lun_vtk)
  
  end subroutine scalin





  !>-------------------------------------------------------------
  !> outputs scalar var in vtk format
  !>
  !> usage:
  !>     call vtkscal(lun,fname,fileformat,time,geo,
  !>                  nvar,utype,uname,uvar)
  !>
  !> where:
  !> \param[in] lun: (integer) unit for error message output
  !> \param[in] fname: (character) filename for main output
  !> \param[in] fileformat: (character) 'ASCII' or 'BINARY'
  !> \param[in] time: (double) simulation time of output
  !> \param[in] geo: (mesh) mesh structure
  !> \param[in] nvar: (integer) number of vars
  !> \param[in] utype: (character array, dim nvar) var types
  !> \param[in] uname: (character array, dim nvar) var names
  !> \param[in] uvar: (double array, dim (nvar,:)) variable
  !>
  !> uses private subroutines alloc_mesh_tmp and
  !>      dealloc_mesh_tmp
  !<-----------------------------------------------------------
  subroutine vtkscal_old(lun,fname,fileformat,time,geo, &
       nvar,utype,uname,uvar)
    use Globals
    use AbstractGeometry
    use LIB_VTK_IO
    implicit none
    ! vars
    integer :: lun
    character(len=*), intent(in) :: fname,fileformat
    real(kind=single), intent(in) :: time
    type(abs_simplex_mesh), intent(in) :: geo
    integer, intent(in) :: nvar
    character(len=*), intent(in) :: utype(nvar),uname(nvar)
    real(kind=double), intent(in) :: uvar(nvar,*)
    ! local vars
    real(kind=double), allocatable :: xcoord(:),ycoord(:),zcoord(:),uuu(:)
    integer, allocatable :: cells(:),ctype(:)    
    character(len=256) :: msg
    character(len=10) :: rdwr
    character(len=4) :: case_arg
    integer :: ncell,nnod,nodes_in_tria
    integer :: res,E_IO
    integer :: i,j
    logical :: var_cell,var_node,rc
    
    ncell = geo%ncell
    nnod = geo%nnode
    nodes_in_tria = geo%nnodeincell
    
    call alloc_mesh_tmp(lun,ncell,nnod,nodes_in_tria,geo,&
         cells,ctype,xcoord,ycoord,zcoord)
    E_IO = VTK_INI(output_format = fileformat, &
         filename = trim(fname), &
         title = 'Unstructured Grid', &
         time = time, &
         mesh_topology = 'UNSTRUCTURED_GRID')
    if (E_IO .ne. 0) rc = IOerr(lun, err_vtk, 'vtkscal', &
         ' error in initializing vtk file, E_IO=',E_IO)
    E_IO = VTK_GEO(Nn = nnod, X = xcoord, Y = ycoord, Z = zcoord)
    if (E_IO .ne. 0) rc = IOerr(lun, err_vtk, 'vtkscal', &
         ' error in writing coordinates in vtk file, E_IO=',E_IO)
    E_IO = VTK_CON(NC = ncell, connect = cells, cell_type = ctype)
    if (E_IO .ne. 0) rc = IOerr(lun, err_vtk, 'vtkscal', &
         ' error in writing mesh topology in vtk file, E_IO=',E_IO)
    
    var_cell=.false.
    var_node=.false.
    do i=1,nvar
       case_arg = trim(adjustl(utype(i)))
       if (case_arg.eq.'cell') var_cell=.true.
       if (case_arg.eq.'node') var_node=.true.
    end do
    if(var_cell) then
       ! write out cell data types (vectors or scalars)
       E_IO = VTK_DAT(NC_NN = ncell, var_location = 'cell')
       if (E_IO .ne. 0) then
          msg = 'writing headers for var ='
          write(rdwr,'(a)') uname(i)
          msg = trim(msg) // trim(rdwr) // ' at component '
          write(rdwr,'(i5)') i
          msg = trim(msg) // trim(rdwr) // ' of '
          write(rdwr,'(i5)') nvar
          msg = trim(msg) // trim(rdwr) // '. E_IO= '
          write(rdwr,'(i5)') E_IO
          msg = trim(msg) // trim(rdwr) 
          rc = IOerr(lun, err_vtk, 'vtkscal', msg, E_IO)
       end if
       allocate(uuu(ncell),stat=res)
       if(res .ne. 0) rc = IOerr(lun, err_alloc, 'vtkscal', &
            'uuu (local array)',res)
       do i=1,nvar
          case_arg = trim(adjustl(utype(i)))
          if(case_arg.eq.'cell') then
             do j=1,ncell
                uuu(j)=uvar(i,j)
             end do
             E_IO = VTK_VAR(NC_NN = ncell, varname = trim(uname(i)), var = uuu )
             if (E_IO .ne. 0) then
                msg = 'writing headers for var ='
                write(rdwr,'(a)') uname(i)
                msg = trim(msg) // trim(rdwr) // ' at component '
                write(rdwr,'(i5)') i
                msg = trim(msg) // trim(rdwr) // ' of '
                write(rdwr,'(i5)') nvar
                msg = trim(msg) // trim(rdwr) // '. E_IO= '
                write(rdwr,'(i5)') E_IO
                msg = trim(msg) // trim(rdwr) 
                rc = IOerr(lun, err_vtk, 'vtkscal', msg, E_IO)
             end if
          end if
       end do
       deallocate(uuu,stat=res)
       if (res.ne.0) rc=IOerr(lun, err_dealloc, 'vtkscal', &
            'uuu (local array)', res)
    end if
    if(var_node) then
       ! write out node data types (vectors or scalars)
       E_IO = VTK_DAT(NC_NN = nnod, var_location = 'node')
       if (E_IO .ne. 0) then
          msg = 'writing headers for var ='
          write(rdwr,'(a)') uname(i)
          msg = trim(msg) // trim(rdwr) // ' at component '
          write(rdwr,'(i5)') i
          msg = trim(msg) // trim(rdwr) // ' of '
          write(rdwr,'(i5)') nvar
          msg = trim(msg) // trim(rdwr) // '. E_IO= '
          write(rdwr,'(i5)') E_IO
          msg = trim(msg) // trim(rdwr) 
          rc = IOerr(lun, err_vtk, 'vtkscal', msg, E_IO)
       end if
       allocate(uuu(nnod),stat=res)
       if(res .ne. 0) rc = IOerr(lun, err_alloc, 'vtkscal', &
            'uuu (local array)',res)
       do i=1,nvar
          case_arg = trim(adjustl(utype(i)))
          if(case_arg.eq.'node') then
             do j=1,nnod
                uuu(j)=uvar(i,j)
             end do
             E_IO = VTK_VAR(NC_NN = nnod, varname = trim(uname(i)), var = uuu )
             if (E_IO .ne. 0) then
                msg = 'writing headers for var ='
                write(rdwr,'(a)') uname(i)
                msg = trim(msg) // trim(rdwr) // ' at component '
                write(rdwr,'(i5)') i
                msg = trim(msg) // trim(rdwr) // ' of '
                write(rdwr,'(i5)') nvar
                msg = trim(msg) // trim(rdwr) // '. E_IO= '
                write(rdwr,'(i5)') E_IO
                msg = trim(msg) // trim(rdwr) 
                rc = IOerr(lun, err_vtk, 'vtkscal', msg, E_IO)
             end if
          end if
       end do
       deallocate(uuu,stat=res)
       if (res.ne.0) rc=IOerr(lun, err_dealloc, 'vtkscal', &
            'uuu (local array)')
    end if
    
    E_IO = VTK_END()
    if (E_IO .ne. 0) then
       write(rdwr,'(i5)') E_IO
       msg = 'error in finalizing vtk file in vtkscal, E_IO='
       msg = trim(msg) // trim(rdwr)
       rc = IOerr(lun, err_vtk, 'vtkscal', msg, E_IO)
    end if
    
    call dealloc_mesh_tmp(lun,cells,ctype,xcoord,ycoord,zcoord) 
    
    return
  end subroutine vtkscal_old


  
  !>-------------------------------------------------------------
  !> allocate temp arrays to be use by vtk subroutines
  !> (private procedure)
  !>
  !> usage:
  !>     call alloc_mesh_tmp(lun,ncell,nnod,nodes_in_tria,geo,
  !>          cells,ctype,xcoord,ycoord,zcoord)
  !>
  !> where:
  !> \param[in] lun: (integer) unit for error message output
  !> \param[in] ncell: (integer) number of cells
  !> \param[in] nnod: (integer) number of nodes
  !> \param[in] nodes_in_tria: (integer) number of nodes per cell
  !> \param[in] geo: (mesh) mesh structure
  !> \param[in] cells: (integer array of dim (ncell*nodes_in_tria)
  !>                   cell topology
  !> \param[in] ctype: (integer array, dim ncell) cell type
  !>               (so far only 5 -> cells)
  !> \param[in] xcoord: (double array, dim nnod) nodal x-coord
  !> \param[in] ycoord: (double array, dim nnod) nodal z-coord
  !> \param[in] zcoord: (double array, dim nnod) nodal z-coord
  !>
  !<-----------------------------------------------------------
  subroutine alloc_mesh_tmp(lun,ncell,nnod,nodes_in_tria,geo,&
       cells,ctype,xcoord,ycoord,zcoord)
    use Globals
    use AbstractGeometry
    use LIB_VTK_IO
    implicit none
    ! vars
    integer :: lun,ncell,nnod,nodes_in_tria
    type(abs_simplex_mesh), intent(in) :: geo
    integer, allocatable, intent(out) :: cells(:),ctype(:)
    real(kind=double), allocatable, intent(out) :: xcoord(:),ycoord(:),zcoord(:)
    ! local vars
    integer :: res
    integer :: i,j,k
    integer :: type_cell
    logical :: rc

    allocate(cells(ncell*(nodes_in_tria+1)),stat=res)
    if(res .ne. 0) rc = IOerr(lun, err_alloc, 'alloc_mesh_tmp', &
         ' cells (local array)',res)
    allocate(ctype(ncell),stat=res)
    if(res .ne. 0) rc = IOerr(lun, err_alloc, 'alloc_mesh_tmp', &
         ' ctype (local array)',res)
    allocate(xcoord(nnod),ycoord(nnod),zcoord(nnod),stat=res)
    if(res .ne. 0) rc = IOerr(lun, err_alloc, 'alloc_mesh_tmp', &
         ' xcoord,ycoord,zcoord (local arrays)',res)

    if ( geo%nnodeincell .eq. 2) type_cell = 3
    if ( geo%nnodeincell .eq. 3) type_cell = 5
    if ( geo%nnodeincell .eq. 4) then
       type_cell = geo%cell_id
    end if

    

    k=0
    do i=1,ncell
       k=k+1
       cells(k) = nodes_in_tria
       do j=1,nodes_in_tria
          k=k+1
          cells(k) = geo%topol(j,i) - 1
       end do
       ctype(i) = type_cell
    end do

    do i=1,nnod
       xcoord(i) = geo%coord(1,i)
       ycoord(i) = geo%coord(2,i)
       zcoord(i) = geo%coord(3,i)
    end do
  end subroutine alloc_mesh_tmp

  !>-------------------------------------------------------------
  !> deallocate temp arrays to be use by vtk subroutines
  !> (private procedure)
  !>
  !> usage:
  !>     call dealloc_mesh_tmp(lun,
  !>          cells,ctype,xcoord,ycoord,zcoord)
  !>
  !> where:
  !> \param[in] lun: (integer) unit for error message output
  !> \param[in] cells: (integer array of dim (ncell*nodes_in_tria)
  !>                   cell topology
  !> \param[in] ctype: (integer array, dim ncell) cell type
  !>               (so far only 5 -> cells)
  !> \param[in] xcoord: (double array, dim nnod) nodal x-coord
  !> \param[in] ycoord: (double array, dim nnod) nodal z-coord
  !> \param[in] zcoord: (double array, dim nnod) nodal z-coord
  !>
  !<-----------------------------------------------------------
  subroutine dealloc_mesh_tmp(lun,cells,ctype,xcoord,ycoord,zcoord) 
    use Globals
    use LIB_VTK_IO
    implicit none
    ! vars
    integer :: lun
    integer, allocatable :: cells(:),ctype(:)
    real(kind=double),allocatable :: xcoord(:),ycoord(:),zcoord(:)
    ! local vars
    integer :: res
    logical :: rc
    
    deallocate(cells,stat=res)
    if (res.ne.0) rc=IOerr(lun, err_dealloc, 'dealloc_mesh_tmp', &
         'cells (local)',res)
    deallocate(ctype,stat=res)
    if (res.ne.0) rc=IOerr(lun, err_dealloc, 'dealloc_mesh_tmp', &
         'ctype (local)',res)
    deallocate(xcoord,ycoord,zcoord,stat=res)
    if (res.ne.0) rc=IOerr(lun, err_dealloc, 'dealloc_mesh_tmp', &
         'xcoord,ycoord,zcoord (local)',res)
  end subroutine dealloc_mesh_tmp

end module vtkloc
