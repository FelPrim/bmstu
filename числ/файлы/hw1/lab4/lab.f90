program lab4
    implicit none
    
    type :: sym_mat
        real, allocatable :: data(:)
        integer :: n = 0
    contains 
        procedure :: init => init_sym
        procedure :: print => print_matrix
        procedure :: get => get_element
        procedure :: set => set_element
    end type

contains
    subroutine init_sym(this, n)
        class(sym_mat), intent(inout) :: this
        integer, intent(in) :: n
        integer :: size

        this%n = n
        size = n*(n+1)/2
        allocate(this%data(size))
        this%data = 0.0
   end subroutine

    function get_element(this, i, j) result(val)
        class(sym_mat), intent(in) :: this
        integer, intent(in) :: this
        real :: val
        integer :: idx

        if (i > j .or. i < 1 .or. j > this%n .or. j < 1) then
            error stop "Indexation error"
        end if
        
        idx = i + j*(j-1)/2
        val = this%data(idx)
    end function

   subroutine set_element(this, i, j, value)
        class(sym_mat), intent(inout) :: this
        integer, intent(in) :: i, j
        real, intent(in) :: value
        integer :: idx

        if (i > j .or. i < 1 .or. j > this%n .or. j < 1) then
            error stop "Indexation error"
        end if

        idx = i + j*(j-1)/2
        this%data(idx) = value
    end subroutine


    subroutine print_matrix(this)
        class(sym_mat), intent(in) :: this
        integer :: i, j

        do i = 1, this % n
            do j = 1, this % n
                if (j < i)
                    write(*, '(F8.2)', advance='no') this%get(j, i)
                else
                    write(*, '(F8.2)', advance='no') this%get(i, j)
                end if
            end do
            print *, ""
        end do
    end subroutine

    subroutine max_nondiag(matrix, alpha, beta, value) 
        class(sym_mat), intent(in) :: this
        integer, intent(out) :: alpha, beta
        real, intent(out) :: value
        integer :: idx
        integer :: size


        size = n*(n+1)/2
        idx = matrix%data(size)
        do i = 1, matrix%n-1
            do


end program 




program simple_matmul
    implicit none
    real, allocatable :: A(:,:), B(:,:), C(:,:)
    integer :: i, j, k, m, n, p
    
    ! Размеры матриц
    m = 3; n = 2; p = 3
    
    ! Выделяем память
    allocate(A(m,n), B(n,p), C(m,p))
    
    ! Инициализация
    A = reshape([1,2,3,4,5,6], [m,n])
    B = reshape([7,8,9,10,11,12], [n,p])
    
    ! Умножение
    do i = 1, m
        do j = 1, p
            C(i,j) = 0
            do k = 1, n
                C(i,j) = C(i,j) + A(i,k) * B(k,j)
            end do
        end do
    end do
    
    ! Вывод
    print *, "Матрица A:"
    do i = 1, m
        print *, A(i,:)
    end do
    
    print *, "Матрица B:" 
    do i = 1, n
        print *, B(i,:)
    end do
    
    print *, "Результат C = A*B:"
    do i = 1, m
        print *, C(i,:)
    end do
    
    deallocate(A, B, C)
end program
