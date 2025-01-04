!!! Codigo para calcular los nodos de Chebyshev y los pesos asociados para la cuadratura de Gauss-Chebyshev
!!! Forma de la integral: f(x) / sqrt(1 - x^2)
!!! https://en.wikipedia.org/wiki/Chebyshev%E2%80%93Gauss_quadrature


!!!!! ESTE MODULO DEBE SER IMPORTADO POR EL PROGRAMA PRINCIPAL !!!!!
module gauss_chebyshev
    implicit none

contains

    subroutine chebyshev_nodes(n, x, w)
        ! Calculo de los nodos y pesos de Chebyshev
        implicit none
        integer, intent(in) :: n
        real(8), intent(out) :: x(n), w(n)
        integer :: i

        do concurrent (i = 1:n) ! Hacer el bucle concurrente para paralelizar todo lo posible 
            x(i) = cos((2.0 * i - 1.0) * 3.141592653589793d0 / (2.0 * n)) ! Nodos de Chebyshev
            w(i) = 3.141592653589793d0 / n ! Peso constante 

        end do
            
    end subroutine chebyshev_nodes

    function gauss(f, x, w, n) result(integral) 
        ! Calculo de la integral con la cuadratura de Gauss-Chebyshev
        implicit none
        real(8), external :: f ! Declarar f como una funcion externa que proviene del programa principal
        integer, intent(in) :: n
        real(8), intent(in) :: x(n), w(n)
        real(8) :: integral
        integer :: i

        integral = 0.0d0 ! Inicializar la integral a 0

        do i = 1, n
            integral = integral + w(i) * f(x(i)) ! Sumar el peso por el valor de la funcion en el nodo
        end do
    end function gauss

end module gauss_chebyshev

