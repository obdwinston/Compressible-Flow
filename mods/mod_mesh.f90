module mod_mesh
    use iso_fortran_env, only: int32, real64
    implicit none

    private
    public :: Mesh

    type :: Element
        real(real64), dimension(2) ::                   coordinates
    end type Element

    type, extends(Element) :: Node
        integer(int32), dimension(20) ::                node_cells
        integer(int32) ::                               n_node_cells
        real(real64), dimension(20) ::                  node_cell_weights
        real(real64), dimension(20) ::                  node_cell_distances
    end type Node

    type, extends(Element) :: Face
        character(20) ::                                face_type
        integer(int32), dimension(2) ::                 face_nodes
        integer(int32), dimension(2) ::                 face_faces ! cell local index
        integer(int32), dimension(2) ::                 face_cells
        real(real64) ::                                 face_area
        real(real64), dimension(2) ::                   face_normal
    end type Face
    
    type, extends(Element) :: Cell
        integer(int32), dimension(3) ::                 cell_nodes
        integer(int32), dimension(3) ::                 cell_faces
        real(real64) ::                                 cell_volume
        real(real64), dimension(3) ::                   cell_face_signs
        integer(int32), dimension(3) ::                 cell_opposite_nodes
        real(real64), dimension(3) ::                   cell_face_distances
        real(real64), dimension(3) ::                   cell_node_distances
    end type Cell

    type :: Mesh
        integer(int32) ::                               n_nodes = 0
        integer(int32) ::                               n_faces = 0
        integer(int32) ::                               n_cells = 0
        type(Node), dimension(:), allocatable ::        nodes
        type(Face), dimension(:), allocatable ::        faces
        type(Cell), dimension(:), allocatable ::        cells
        character(20), dimension(:), allocatable ::     face_types
        integer(int32), dimension(:), allocatable ::    n_face_types
        integer(int32), dimension(:), allocatable ::    surface_nodes
    contains
        procedure, pass(self) :: set_mesh
        procedure, pass(self) :: save_mesh
        procedure, pass(self) :: load_mesh
    end type Mesh
    
contains

    subroutine set_mesh(self, mesh_file, body_file)
        class(Mesh), intent(in out) :: self
        character(*), intent(in) :: mesh_file, body_file
        open(10, file=trim(mesh_file), action='read')
        call set_cells(self, 10)
        call set_nodes(self, 10)
        call set_faces(self, 10)
        close(10)
        call set_surface_nodes(self, body_file)
        call set_cell_faces(self)
        call set_face_cells(self)
        call set_cell_components(self) ! order important
        call set_face_components(self) ! order important
        call set_node_components(self) ! order important
        call set_gradient_components(self)
    end subroutine set_mesh

    subroutine save_mesh(self, save_file)
        class(Mesh), intent(in) :: self
        character(*), intent(in) :: save_file
        open(10, file=trim(save_file), action='write')
        write(10, *) self % n_nodes, self % n_faces, self % n_cells
        write(10, *) self % nodes
        write(10, *) self % faces
        write(10, *) self % cells
        write(10, *) size(self % face_types)
        write(10, *) self % face_types
        write(10, *) self % n_face_types
        write(10, *) size(self % surface_nodes)
        write(10, *) self % surface_nodes
        close(10)
    end subroutine save_mesh
    
    subroutine load_mesh(self, load_file)
        class(Mesh), intent(in out) :: self
        character(*), intent(in) :: load_file
        integer(int32) :: k
        open(10, file=trim(load_file), action='read')
        read(10, *) self % n_nodes, self % n_faces, self % n_cells
        allocate(self % nodes(self % n_nodes))
        allocate(self % faces(self % n_faces))
        allocate(self % cells(self % n_cells))
        read(10, *) self % nodes
        read(10, *) self % faces
        read(10, *) self % cells
        read(10, *) k
        allocate(self % face_types(k))
        allocate(self % n_face_types(k))
        read(10, *) self % face_types
        read(10, *) self % n_face_types
        read(10, *) k
        allocate(self % surface_nodes(k))
        read(10, *) self % surface_nodes
        close(10)
    end subroutine load_mesh

    subroutine set_cells(self, file_unit)
        type(Mesh), intent(in out) :: self
        integer(int32), intent(in) :: file_unit
        character(100) :: line
        integer(int32) :: i, n1, n2, n3
        read(10, *)
        read(10, *) line, self % n_cells
        allocate(self % cells(self % n_cells))
        do i = 1, self % n_cells
            read(10, *) line, n1, n2, n3
            self % cells(i) % cell_nodes(1) = n1 + 1
            self % cells(i) % cell_nodes(2) = n2 + 1
            self % cells(i) % cell_nodes(3) = n3 + 1
        end do
    end subroutine set_cells

    subroutine set_nodes(self, file_unit)
        type(Mesh), intent(in out) :: self
        integer(int32), intent(in) :: file_unit
        character(100) :: line
        integer(int32) :: i
        real(real64) :: nx, ny
        read(10, *) line, self % n_nodes
        allocate(self % nodes(self % n_nodes))
        do i = 1, self % n_nodes
            read(10, *) nx, ny
            self % nodes(i) % coordinates = [nx, ny]
        end do
    end subroutine set_nodes

    subroutine set_faces(self, file_unit)
        type(Mesh), intent(in out) :: self
        integer(int32), intent(in) :: file_unit
        character(100) :: line
        integer(int32) :: i, k
        character(20), dimension(:), allocatable :: c1
        integer(int32), dimension(:), allocatable :: l1, l2
        c1 = [character(20) ::] ! type list
        l1 = [integer(int32) ::] ! node list
        l2 = [integer(int32) ::] ! node list
        read(10, *) line, k
        allocate(self % face_types(k + 1))
        allocate(self % n_face_types(k + 1))
        call set_boundary_faces(self, c1, l1, l2, file_unit, k)
        call set_interior_faces(self, c1, l1, l2, file_unit)
        self % n_faces = size(l1)
        allocate(self % faces(size(l1)))
        do i = 1, self % n_faces
            self % faces(i) % face_type = c1(i)
            self % faces(i) % face_nodes(1) = l1(i)
            self % faces(i) % face_nodes(2) = l2(i)
        end do
    end subroutine set_faces

    subroutine set_boundary_faces(self, c1, l1, l2, file_unit, k)
        type(Mesh), intent(in out) :: self
        character(20), dimension(:), allocatable, intent(in out) :: c1
        integer(int32), dimension(:), allocatable, intent(in out) :: l1, l2
        integer(int32), intent(in) :: file_unit, k
        character(100) :: line
        integer(int32) :: i, j, n1, n2
        do i = 1, k
            read(10, *) line, self % face_types(i)
            read(10, *) line, self % n_face_types(i)
            do j = 1, self % n_face_types(i)
                read(10, *) line, n1, n2
                c1 = [c1, self % face_types(i)]
                l1 = [l1, n1 + 1]
                l2 = [l2, n2 + 1]
            end do
        end do
    end subroutine set_boundary_faces

    subroutine set_interior_faces(self, c1, l1, l2, file_unit)
        type(Mesh), intent(in out) :: self
        character(20), dimension(:), allocatable, intent(in out) :: c1
        integer(int32), dimension(:), allocatable, intent(in out) :: l1, l2
        integer(int32), intent(in) :: file_unit
        integer(int32) :: i, n1, n2, n3
        do i = 1, self % n_cells
            n1 = self % cells(i) % cell_nodes(1)
            n2 = self % cells(i) % cell_nodes(2)
            n3 = self % cells(i) % cell_nodes(3)
            call set_face(c1, l1, l2, n1, n2)
            call set_face(c1, l1, l2, n2, n3)
            call set_face(c1, l1, l2, n3, n1)
        end do
        self % face_types(size(self % face_types)) = 'INTERIOR'
        self % n_face_types(size(self % n_face_types)) = &
        size(l1) - sum(self % n_face_types)
    end subroutine set_interior_faces

    subroutine set_surface_nodes(self, body_file)
        type(Mesh), intent(in out) :: self
        character(*), intent(in) :: body_file
        integer(int32) :: i, ios, n
        integer(int32), dimension(:), allocatable :: nodes
        real(real64) :: nx, ny
        nodes = [integer(int32) ::]
        open(11, file=trim(body_file), action='read', iostat=ios)
        do
            read(11, *, iostat=ios) nx, ny
            if (ios /= 0) exit
            n = get_node(self, nx, ny)
            nodes = [nodes, n]
        end do
        allocate(self % surface_nodes, source=nodes)
        close(11)
    end subroutine set_surface_nodes

    subroutine set_cell_faces(self)
        type(Mesh), intent(in out) :: self
        integer(int32) :: i, n1, n2, n3
        integer(int32), dimension(:), allocatable :: l1, l2
        allocate(l1(self % n_faces))
        allocate(l2(self % n_faces))
        l1 = self % faces(:) % face_nodes(1)
        l2 = self % faces(:) % face_nodes(2)
        do i = 1, self % n_cells
            n1 = self % cells(i) % cell_nodes(1)
            n2 = self % cells(i) % cell_nodes(2)
            n3 = self % cells(i) % cell_nodes(3)
            self % cells(i) % cell_faces(1) = get_face(l1, l2, n1, n2)
            self % cells(i) % cell_faces(2) = get_face(l1, l2, n2, n3)
            self % cells(i) % cell_faces(3) = get_face(l1, l2, n3, n1)
        end do
    end subroutine set_cell_faces

    subroutine set_face_cells(self)
        type(Mesh), intent(in out) :: self
        integer(int32) :: i, j
        integer(int32), dimension(:), allocatable :: l3
        do i = 1, self % n_faces
            l3 = [integer(int32) ::] ! cell list
            do j = 1, self % n_cells
                if (any(self % cells(j) % cell_faces == i)) l3 = [l3, j]
                if (size(l3) == 2) exit
            end do
            if (size(l3) == 2) then
                self % faces(i) % face_cells(1) = l3(1)
                self % faces(i) % face_cells(2) = l3(2)
            else
                self % faces(i) % face_cells(1) = l3(1)
                self % faces(i) % face_cells(2) = l3(1)
            end if
        end do
    end subroutine set_face_cells

    subroutine set_face_components(self)
        type(Mesh), intent(in out) :: self
        integer(int32) :: i, j, cj, n1, n2
        real(real64) :: n1x, n1y, n2x, n2y
        do i = 1, self % n_faces
            n1 = self % faces(i) % face_nodes(1)
            n2 = self % faces(i) % face_nodes(2)
            n1x = self % nodes(n1) % coordinates(1)
            n1y = self % nodes(n1) % coordinates(2)
            n2x = self % nodes(n2) % coordinates(1)
            n2y = self % nodes(n2) % coordinates(2)
            self % faces(i) % coordinates = &
            [0.5*(n2x + n1x), 0.5*(n2y + n1y)]
            self % faces(i) % face_area = &
            sqrt((n2x - n1x)**2 + (n2y - n1y)**2)
            self % faces(i) % face_normal = &
            [(n2y - n1y)/self % faces(i) % face_area, &
            -(n2x - n1x)/self % faces(i) % face_area]
            do j = 1, 2
                cj = self % faces(i) % face_cells(j)
                self % faces(i) % face_faces(j) = &
                findloc(self % cells(cj) % cell_faces, i, 1)
            end do
        end do
    end subroutine set_face_components

    subroutine set_cell_components(self)
        type(Mesh), intent(in out) :: self
        integer(int32) :: i, j, n1, n2, n3, fj, cj
        real(real64) :: n1x, n1y, n2x, n2y, n3x, n3y
        do i = 1, self % n_cells
            n1 = self % cells(i) % cell_nodes(1)
            n2 = self % cells(i) % cell_nodes(2)
            n3 = self % cells(i) % cell_nodes(3)
            n1x = self % nodes(n1) % coordinates(1)
            n1y = self % nodes(n1) % coordinates(2)
            n2x = self % nodes(n2) % coordinates(1)
            n2y = self % nodes(n2) % coordinates(2)
            n3x = self % nodes(n3) % coordinates(1)
            n3y = self % nodes(n3) % coordinates(2)
            self % cells(i) % coordinates = &
            [(n1x + n2x + n3x)/3, (n1y + n2y + n3y)/3]
            self % cells(i) % cell_volume = &
            0.5*abs(n1x*(n2y - n3y) + n2x*(n3y - n1y) + n3x*(n1y - n2y))
            do j = 1, 3
                fj = self % cells(i) % cell_faces(j)
                if (i == self % faces(fj) % face_cells(1)) then
                    self % cells(i) % cell_face_signs(j) = 1.0
                else
                    self % cells(i) % cell_face_signs(j) = -1.0
                end if
            end do
        end do
    end subroutine set_cell_components

    subroutine set_gradient_components(self)
        type(Mesh), intent(in out) :: self
        integer(int32) :: i, j, n1, n2, n3, fj, na, nb, nc
        real(real64) :: cx, cy, fx, fy, nx, ny
        do i = 1, self % n_cells
            n1 = self % cells(i) % cell_nodes(1)
            n2 = self % cells(i) % cell_nodes(2)
            n3 = self % cells(i) % cell_nodes(3)
            do j = 1, 3
                fj = self % cells(i) % cell_faces(j)
                na = self % faces(fj) % face_nodes(1)
                nb = self % faces(fj) % face_nodes(2)
                if (n1 /= na .and. n1 /= nb) nc = n1
                if (n2 /= na .and. n2 /= nb) nc = n2
                if (n3 /= na .and. n3 /= nb) nc = n3
                self % cells(i) % cell_opposite_nodes(j) = nc
                cx = self % cells(i) % coordinates(1)
                cy = self % cells(i) % coordinates(2)
                fx = self % faces(fj) % coordinates(1)
                fy = self % faces(fj) % coordinates(2)
                nx = self % nodes(nc) % coordinates(1)
                ny = self % nodes(nc) % coordinates(2)
                self % cells(i) % cell_face_distances(j) = &
                sqrt((fx - cx)**2 + (fy - cy)**2)
                self % cells(i) % cell_node_distances(j) = &
                sqrt((cx - nx)**2 + (cy - ny)**2)
            end do
        end do
    end subroutine set_gradient_components

    subroutine set_node_components(self)
        type(Mesh), intent(in out) :: self
        integer(int32) :: i, j, k
        real(real64) :: nx, ny, cx, cy
        do i = 1, self % n_nodes
            k = 1 ! must start from 1
            do j = 1, self % n_cells
                if (any(self % cells(j) % cell_nodes == i)) then
                    nx = self % nodes(i) % coordinates(1)
                    ny = self % nodes(i) % coordinates(2)
                    cx = self % cells(j) % coordinates(1)
                    cy = self % cells(j) % coordinates(2)
                    self % nodes(i) % node_cells(k) = j
                    self % nodes(i) % node_cell_distances(k) = &
                    1/sqrt((nx - cx)**2 + (ny - cy)**2)
                    k = k + 1
                end if
            end do
            self % nodes(i) % n_node_cells = k - 1
            do j = 1, self % nodes(i) % n_node_cells
                self % nodes(i) % node_cell_weights(j) = &
                self % nodes(i) % node_cell_distances(j)/ &
                sum(self % nodes(i) % node_cell_distances)
            end do
        end do
    end subroutine set_node_components

    subroutine set_face(c1, l1, l2, n1, n2)
        character(20), dimension(:), allocatable, intent(in out) :: c1
        integer(int32), dimension(:), allocatable, intent(in out) :: l1, l2
        integer(int32), intent(in) :: n1, n2
        integer(int32), dimension(:), allocatable :: m1, m2
        logical :: in
        allocate(m1, mold=l1)
        allocate(m2, mold=l2)
        m1 = n1 ! array of n1 with length l1
        m2 = n2 ! array of n2 with length l2
        in = any((l1 == m1) .and. (l2 == m2)) .or. any((l1 == m2) .and. (l2 == m1))
        if (.not. in) then
            c1 = [c1, 'INTERIOR']
            l1 = [l1, n1]
            l2 = [l2, n2]
        end if
    end subroutine set_face

    pure function get_face(l1, l2, n1, n2) result(res)
        integer(int32), dimension(:), allocatable, intent(in) :: l1, l2
        integer(int32), intent(in) :: n1, n2
        integer(int32) :: res
        do res = 1, size(l1)
            if ((l1(res) == n1 .and. l2(res) == n2) .or. &
            (l1(res) == n2 .and. l2(res) == n1)) return
        end do
    end function get_face

    function get_node(msh, nx, ny) result(res)
        type(Mesh), intent(in) :: msh
        real(real64), intent(in) :: nx, ny
        integer(int32) :: res
        real(real64) :: mx, my, dx, dy
        do res = 1, msh % n_nodes
            mx = msh % nodes(res) % coordinates(1)
            my = msh % nodes(res) % coordinates(2)
            dx = abs(mx - (nx + 24.5))
            dy = abs(my - (ny + 25.))
            if (dx <= 1e-5 .and. dy <= 1e-5) return
        end do
    end function get_node
    
end module mod_mesh
