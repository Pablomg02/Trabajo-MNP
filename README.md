# Guia de uso
## Introducción
Este código permite llevar a cabo la optimización de la distribución de sustentación en un ala.
- En la carpeta `src` se encuentra todo el codigo fuente utilizado
- En la carpeta `kernel` se encuentra el codigo Fortran necesario para llevar a cabo la integración.

## Instalación
Dado que se trata de una libreria de Fortran en Python y no cuenta con un instalador de la libreria automatico, es necesario
seguir los siguientes pasos en la terminal para utilizarla.

1. **Clonar la libreria** (en caso de estar trabajando desde local o en una plataforma que no es CoCalc):
```bash
git clone https://github.com/Pablomg02/Trabajo-MNP.git
```

2. Dirigirse a la ruta donde se encuentra el archivo fortran (este comando puede cambiar en funcion de la ruta inicial
```bash
cd src/kernel
```

3. Una vez en esa dirección, es necesario tener la librería Numpy instalada, ya que esta trae incluido f2py, una herramienta
que permite compilar archivos de Fortran para Python fácilmente. Normalmente, Numpy ya se encuentra instalado, pero si hiciese
falta, se instala con el siguiente comando
```bash
pip install numpy
```

4. Finalmente, se compila el archivo con el siguiente comando (**Es importante estar en la direccion correcta y tener NumPy instalado**)
```bash
f2py -c -m gauss_chebyshev gauss_chebyshev.f90
```

5. Comprobar que el archivo compilado está en la misma carpeta que `gauss_chebyshev.f90` y tiene el mismo nombre

Una vez realizado esto, el código en `main.ipynb` es completamente funcional
