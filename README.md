# Electrodinámica 2026: Método de los Momentos (MoM) y Elementos Finitos (MFEM)

Este repositorio contiene la implementación numérica para la resolución de problemas electrostáticos tridimensionales mediante el **Método de los Momentos (MoM)** en Python y el **Método de Elementos Finitos (FEM)** en C++ utilizando la biblioteca **MFEM**.


---

## 📁 Estructura del Proyecto

```text
Electrodynamics2026/
│
├── README.md                      # Documentación y guía de ejecución
├── requirements.txt               # Lista de dependencias de Python
├── main.py                        # Script principal ejecutable para generar todas las gráficas
│
├── src/                           # Módulo fuente principal en Python (Estructura limpia)
│   ├── __init__.py                # Exportación del paquete
│   ├── green_functions.py         # Funciones de Green discreta (Point-Matching y Galerkin)
│   ├── single_plate.py            # Solucionador MoM para placa conductora aislada
│   ├── dual_plates.py             # Solucionador MoM para dos placas paralelas de tamaño arbitrario
│   ├── cylinder.py                # Solucionador MoM para cilindro conductor cerrado ante cargas
│   └── utils.py                   # Utilidades de gestión de directorio output/ y figuras
│
├── output/                        # Carpeta donde se guardan automáticamente las gráficas (.png)
│   ├── 01_green_functions_coplanar.png
│   ├── 02_green_functions_parallel.png
│   ├── 03_green_functions_surface_h03.png
│   ├── 04_single_plate_charge_density.png
│   ├── 05_dual_plates_charge_density.png
│   ├── 06_cylinder_charge_density.png
│   └── 07_cylinder_potential_contours.png
│
├── notebooks/                     # Cuadernos Jupyter para análisis interactivo
│   └── placas_distintos_tamanos.ipynb
│
└── mfem/                          # Código C++ con la biblioteca MFEM
    ├── docker-compose.yml         # Contenedor Docker listo para compilar/ejecutar MFEM
    ├── poisson_constant/          # Ecuación de Poisson con fuente constante
    │   ├── poisson.cpp
    │   └── makefile
    └── poisson_delta/             # Ecuación de Poisson con fuente puntual (Delta de Dirac)
        ├── Makefile
        ├── apps/
        │   └── main.cpp
        ├── include/
        └── src/
```

---

## Guía de ejecución

Para reproducir todos los análisis y generar automáticamente las 7 gráficas principales en la carpeta `output/`, siga los siguientes pasos:

### 1. Instalación de Dependencias
Instale los paquetes requeridos (`numpy`, `matplotlib`, `scipy`) ejecutando en terminal:

```bash
pip install -r requirements.txt
```

### 3. Ejecución del Análisis Principal
Ejecute el script principal `main.py`:

```bash
python main.py
```

Al finalizar la ejecución, se mostra en consola los cálculos numéricos de capacitancia y reciprocidad, y todas las figuras resultantes quedan guardadas automáticamente en la carpeta **`output/`**.

---

## 📊 Gráficas Generadas en `output/`

| Archivo | Descripción del Análisis Electroestático |
| :--- | :--- |
| `01_green_functions_coplanar.png` | Comparación de las Funciones de Green discretas para Galerkin ($\Psi$), Point-Matching ($\Gamma$) y Espacio Libre ($G$) en celdas coplanares ($h=0$). |
| `02_green_functions_parallel.png` | Comparación de las Funciones de Green para celdas paralelas variando la distancia vertical $z_c = h$. |
| `03_green_functions_surface_h03.png` | Superficie 3D de la interacción entre celdas a una separación $h = 0.3$. |
| `04_single_plate_charge_density.png` | Distribución 3D de la densidad de carga inducida $\rho_s(x,y)$ en una placa aislada a potencial constante. |
| `05_dual_plates_charge_density.png` | Densidad de carga 3D en un sistema de dos placas de distinto tamaño ($A_{\text{top}}=1.0$, $A_{\text{bot}}=0.5$) y verificación de reciprocidad ($C_{12} = C_{21}$). |
| `06_cylinder_charge_density.png` | Densidad de carga inducida en la superficie de un cilindro conductor cerrado sometido a cargas puntuales interna y externa. |
| `07_cylinder_potential_contours.png` | Mapas de contorno 2D del potencial total $V(x,y)$ evaluados en el plano de observación $z = 0.5$. |

---

## 💻 Módulo C++ MFEM (Elementos Finitos) (opcional)

La carpeta `mfem/` contiene dos aplicaciones en C++ diseñadas para resolver la Ecuación de Poisson $\nabla^2 V = -f$:

1. **`mfem/poisson_constant/`**: Resuelve la ecuación con un término de fuente constante $f(x,y)=1$.
2. **`mfem/poisson_delta/`**: Resuelve la ecuación con una fuente puntual (Delta de Dirac).

### Compilación y Ejecución con Docker
Si dispone de Docker y Docker Compose, puede compilar e interactuar con MFEM en un entorno listo para usar:

```bash
cd mfem
docker-compose up -d
docker exec -it mfem-dev bash
```

Dentro del contenedor, navegue a `poisson_constant` o `poisson_delta` y ejecute `make` para compilar la aplicación C++.

---