!!! Codigo para calcular los nodos de Chebyshev y los pesos asociados para la cuadratura de Gauss-Chebyshev
!!! Forma de la integral: f(x) / sqrt(1 - x^2)
!!! https://en.wikipedia.org/wiki/Chebyshev%E2%80%93Gauss_quadrature


    function run(n) result(result)
        
        implicit none
        integer, intent(in) :: n
        
        real(8) :: result
        real(8) :: x(n), w(n)
        real(8) :: gauss

        call chebyshev_nodes(n, x, w)
        result = gauss(x, w, n)

    end function run



    subroutine chebyshev_nodes(n, x, w)
        ! Calculo de los nodos y pesos de Chebyshev
        implicit none
        integer, intent(in) :: n
        real(8), intent(out) :: x(n), w(n)
        real(8), parameter :: pi = 3.141592653589793d0
        integer :: i

        do concurrent (i = 1:n)
            x(i) = cos(real(i, 8) / real(n + 1, 8) * pi)
            w(i) = (pi / real(n + 1, 8)) * sin(real(i, 8) / real(n + 1, 8) * pi)**2
        end do
    end subroutine chebyshev_nodes




    function gauss(x, w, n) result(integral)
        ! Calculo de la integral con la cuadratura de Gauss-Chebyshev
        implicit none
        integer, intent(in) :: n
        real(8), intent(in) :: x(n), w(n)
        real(8) :: integral
        integer :: i
        real(8) :: f

        integral = 0.0d0

        do i = 1, n
            integral = integral + w(i) * f(x(i))
        end do
    end function gauss



    function f(x) result(result)
        ! Definición de la función polinómica f(x) con coeficientes en a
        implicit none
        real(8), intent(in) :: x
        real(8) :: result
        
        result = x**2
    end function f

