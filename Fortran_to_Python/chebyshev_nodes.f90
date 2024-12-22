!!! Codigo para calcular los nodos de Chebyshev y los pesos asociados para la cuadratura de Gauss-Chebyshev
!!! Forma de la integral: f(x) / sqrt(1 - x^2)
!!! https://en.wikipedia.org/wiki/Chebyshev%E2%80%93Gauss_quadrature

    subroutine chebyshev_nodes(n, x, w)
        ! Calculo de los nodos y pesos de Chebyshev
        implicit none
        integer, intent(in) :: n
        real(8), intent(out) :: x(n), w(n)
        integer :: i

        do concurrent (i = 1:n)
            !x(i) = 0.5 * (a + b) + 0.5 * (b - a) * cos((2.0 * i - 1.0) * 3.141592653589793d0 / (2.0 * n))
            x(i) = cos((2.0 * i - 1.0) * 3.141592653589793d0 / (2.0 * n))
            w(i) = 3.141592653589793d0 / n ! Peso constante

            !print *, "x(", i, ") = ", x(i), " ; w(", i, ") = ", w(i)
        end do
            
    end subroutine chebyshev_nodes

