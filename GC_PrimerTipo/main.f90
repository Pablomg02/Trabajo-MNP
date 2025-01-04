!!!!!!! PROGRAMA PRINCIPAL !!!!!!!
!!!!!!! HACE USO DEL MODULO GAUSS_CHEBYSHEV !!!!!!!


program main
    use gauss_chebyshev ! Importar el modulo con las subrutinas
    implicit none

    integer :: n, i ! Declarar variables enteras
    real(8), allocatable :: x(:), w(:) ! Declarar arreglos de reales, que se asignaran en tiempo de ejecucion
    real(8) :: result 
    real :: start_time, end_time

    print *, "Introduce el numero de nodos n:"
    read *, n

    call cpu_time(start_time) ! Iniciar el contador de tiempo

    allocate(x(n)) ! Asignar memoria para el arreglo x, de tamaño ahora sabido
    allocate(w(n))

    ! Llamar a la subrutina para calcular los nodos y pesos
    call chebyshev_nodes(n, x, w)
    
    ! Calcular la integral con el modulo gauss_chebyshev
    result = gauss(f, x, w, n) ! f es la funcion definida mas abajo, y aparece en f(x)/sqrt(1-x^2)

    call cpu_time(end_time) ! Finalizar el contador de tiempo

    ! Liberar memoria
    deallocate(x)
    deallocate(w)

    print *, "Tiempo de calculo: ", end_time - start_time 

    print *, "El resultado de la integral es: ", result


    contains 
        function f(x) result(result) ! Aqui se define la funcion f(x) que se va a integrar
            implicit none
            real(8), intent(in) :: x
            real(8) :: result
        
            ! Cálculo de la función
            result = cos(x)
        end function f


end program main


