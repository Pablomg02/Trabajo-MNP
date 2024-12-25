!!! Codigo para calcular los nodos de Chebyshev y los pesos asociados para la cuadratura de Gauss-Chebyshev
!!! Forma de la integral: f(x) * sqrt(1 - x^2)
!!! https://en.wikipedia.org/wiki/Chebyshev%E2%80%93Gauss_quadrature


    module gauss_chebyshev
    implicit none

    contains

        function run_lift(n, a) result(res)
            implicit none
            integer, intent(in) :: n
            real(8) :: res
            real(8) :: x(n), w(n)
            real(8), intent(in) :: a(3)

            call chebyshev_nodes(n, x, w)
            res = gauss(x, w, n, a)
        end function run_lift

        function run_mom(n, a) result(res)
            implicit none
            integer, intent(in) :: n
            real(8) :: res
            real(8) :: x(n), w(n)
            real(8), intent(in) :: a(3)

            call chebyshev_nodes(n, x, w)
            res = gauss_x(x, w, n, a) ! Usa la funcion multiplicada por x
        end function run_mom


        
        subroutine chebyshev_nodes(n, x, w)
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



        function gauss(x, w, n, a) result(integral)
            implicit none
            integer, intent(in) :: n
            real(8), intent(in) :: x(n), w(n)
            real(8), intent(in) :: a(3)
            real(8) :: integral
            integer :: i

            integral = 0.0d0
            do i = 1, n
                integral = integral + w(i) * f(x(i), a)
            end do
        end function gauss


        function gauss_x(x, w, n, a) result(integral)
            implicit none
            integer, intent(in) :: n
            real(8), intent(in) :: x(n), w(n)
            real(8), intent(in) :: a(3)
            real(8) :: integral
            integer :: i

            integral = 0.0d0
            do i = 1, n
                integral = integral + w(i) * f_x(x(i), a)
            end do
        end function gauss_x



        function f(x, a) result(res)
            implicit none
            real(8), intent(in) :: x
            real(8), intent(in) :: a(3)
            real(8) :: res

            res = a(1) * x + a(2) * x**2 + a(3) * x**3

        end function f

        function f_x(x, a) result(res)
            implicit none
            real(8), intent(in) :: x
            real(8), intent(in) :: a(3)
            real(8) :: res

            res = (a(1) * x + a(2) * x**2 + a(3) * x**3) * (x+1)

        end function f_x




    end module gauss_chebyshev

