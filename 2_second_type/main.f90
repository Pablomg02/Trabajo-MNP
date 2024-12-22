program main
    use gauss_chebyshev
    implicit none

    integer :: n, i
    real(8), allocatable :: x(:), w(:)
    real(8) :: result
    real :: start_time, end_time

    ! Solicitar los valores de a, b y n al usuario
    !print *, "Introduce el valor de a:"
    !read *, a
    !print *, "Introduce el valor de b:"
    !read *, b
    print *, "Introduce el numero de nodos n:"
    read *, n

    call cpu_time(start_time)

    allocate(x(n))
    allocate(w(n))

    ! Llamar a la subrutina para calcular los nodos y pesos
    call chebyshev_nodes(n, x, w)
    
    ! Calcular la integral
    result = gauss(f, x, w, n)

    call cpu_time(end_time)

    ! Liberar memoria
    deallocate(x)
    deallocate(w)

    print *, "Tiempo de calculo: ", end_time - start_time

    print *, "El resultado de la integral es: ", result


    contains
        function f(x) result(result)
            implicit none
            real(8), intent(in) :: x
            real(8) :: result
        
            ! Cálculo de la función
            result = 3*x**2 + 2*x
        end function f


end program main


