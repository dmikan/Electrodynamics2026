import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# --- 1. FUNCIONES FUNDAMENTALES ---
def definite(prim, x1, x2, y1, y2):
    """
    Calcula la integral doble de una función f(x,y) sobre un dominio 
    rectangular definido por [x1, x2] y [y1, y2].
    
    Esta función implementa una generalización del Teorema Fundamental del
    Cálculo para dos dimensiones.
    
    El valor de la integral se obtiene evaluando dicha primitiva en las cuatro 
    esquinas del rectángulo.
    """
    return prim(x1, y1) + prim(x2, y2) - prim(x1, y2) - prim(x2, y1)

def DiscrGF(xc, yc, zc, a, b):
    """
    Calcula la Función de Green Discreta para Point-Matching.
    """
    def prim1(u, v, h):
        return ((u * np.arcsinh(v/np.sqrt(u**2 + h**2))) +
                (v * np.arcsinh(u/np.sqrt(v**2 + h**2))) -
                (h * np.arctan(u*v/(h*np.sqrt(u**2 + v**2 + h**2))))
                )
    # Pequeño valor para evitar divisiones por cero en el plano z=0
    zc = zc + np.finfo(float).eps 
    """
    Evaluación de la primitiva en las 4 esquinas de la celda fuente.
    
    Como la celda de fuente está centrada en el origen, los límites de
    integración van de −a/2 a +a/2 y de −b/2 a +b/2.
    """
    # Retorna el potencial normalizado por el área de la celda
    result = (prim1(xc + a/2, yc + b/2, zc) + prim1(xc - a/2, yc - b/2, zc) -
             prim1(xc + a/2, yc - b/2, zc) - prim1(xc - a/2, yc + b/2, zc))
    return result/(a*b)

def AvDiscrGF(xc, yc, zc, a, b):
    """
    Calcula la interacción promedio entre dos celdas rectangulares (Galerkin).
    
    Implementa la reducción de la integral 4D a 2D usando 4 primitivas 
    analíticas.
    """
    # Evitar división por cero
    zc = zc + np.finfo(float).eps  
    # Definición de las 4 primitivas de la Tabla 1 del artículo
    prim1 = lambda u, v: ((u * np.arcsinh(v/np.sqrt(u**2 + zc**2))) +
                          (v * np.arcsinh(u/np.sqrt(v**2 + zc**2))) -
                          (zc * np.arctan(u*v/(zc*np.sqrt(u**2 + v**2 + zc**2))))
                          )
    prim2 = lambda u, v: (1/2*(v*np.sqrt(u**2 + v**2 + zc**2)
                         + (u**2 + zc**2)*np.arcsinh(v/np.sqrt(u**2 + zc**2))))
    prim3 = lambda u, v: (1/2*(u*np.sqrt(u**2 + v**2 + zc**2)
                         + (v**2 + zc**2)*np.arcisnh(u/np.sqrt(v**2 + zc**2))))
    prim4 = lambda u, v: (1/3*(u**2 + v**2 + zc**2)**(3/2))
    
    # Cálculo de las contribuciones en los diferentes límites (Yuu, Yud, Ydu, Ydd)
    Yuu = ((xc - a)*(yc - b)*definite(prim1, xc-a, xc, yc-b, yc)
           - (yc - b)*definite(prim2, xc-a, xc, yc-b, yc)
           - (xc - a)*definite(prim3, xc-a, xc, yc-b, yc)
           + definite(prim4, xc-a, xc, yc-b, yc))
    Yud = ((xc - a)*(yc + b)*definite(prim1, xc-a, xc, yc, yc+b)
           + (yc + b)*definite(prim2, xc-a, xc, yc, yc+b)
           + (xc - a)*definite(prim3, xc-a, xc, yc, yc+b)
           - definite(prim4, xc-a, xc, yc, yc+b))
    Ydu = ((xc + a)*(yc - b)*definite(prim1, xc, xc+a, yc-b, yc)
           + (yc - b)*definite(prim2, xc, xc+a, yc-b, yc)
           + (xc + a)*definite(prim3, xc, xc+a, yc-b, yc)
           - definite(prim4, xc, xc+a, yc-b, yc))
    Ydd = ((xc + a)*(yc + b)*definite(prim1, xc, xc+a, yc, yc+b)
           - (yc + b)*definite(prim2, xc, xc+a, yc, yc+b)
           - (xc + a)*definite(prim3, xc, xc+a, yc, yc+b)
           + definite(prim4, xc, xc+a, yc, yc+b))
    # El resultado final se normaliza por el producto de las áreas de ambas celdas
    return (Yuu + Yud + Ydu + Ydd) / (a * b)**2


# --- 2. SOLUCION PLACA ESTÁTICA ---
def solve_plate_vectorized():
    # Configuración de la placa [1]
    aside, bside = 1.0, 1.0
    M, N = 20, 20  # Divisiones por lado
    a, b = aside / M, bside / N
    tresh = 10     # Umbral de aproximación
    flag_g = 0     # 0 para Point-Matching, 1 para Galerkin

    # Mallado de centros de celdas
    x = np.linspace(a/2, aside - a/2, M)
    y = np.linspace(b/2, bside - b/2, N)
    xmat, ymat = np.meshgrid(x, y)
    
    # Aplanado de vectores (order='F' para coincidir con Matlab)
    xvec = xmat.flatten(order='F')
    yvec = ymat.flatten(order='F')
    ndim = M * N

    # --- VECTORIZACIÓN: Llenado de la matriz MoM ---
    # xvec[:, None] crea un vector columna, xvec[None, :] crea uno fila
    XC = np.abs(xvec[:, np.newaxis] - xvec[np.newaxis, :])
    YC = np.abs(yvec[:, np.newaxis] - yvec[np.newaxis, :])
    ZC = 0.0 
    
    R = np.sqrt(XC**2 + YC**2 + ZC**2)
    rh = np.sqrt(XC**2 + YC**2)

    # Selección de estrategia basada en el umbral
    # Si rh > tresh * max(a,b), usamos 1/R. Si no, usamos la función discreta.
    threshold_val = tresh * max(a, b)
    
    if flag_g == 1:
        mom = np.where(rh > threshold_val, 1/R, AvDiscrGF(XC, YC, ZC, a, b))
    else:
        mom = np.where(rh > threshold_val, 1/R, DiscrGF(XC, YC, ZC, a, b))

    # 3) Vector de excitacion
    
    # ESCENARIO 1: Potencial constante unitario
    unorm = 1.0
    #excv = unorm * np.ones((ndim, 1))
    
    # ESCENARIO 2: Carga puntual unitaria
    xq, yq, zq = 0.5, 0.5, 0.5  # Posición de la carga puntual sobre el centro
    excv = 1.0 / np.sqrt((xq - xvec)**2 + (yq - yvec)**2 + zq**2)
    excv = excv.reshape(-1, 1)

    # 4) Resolución y cálculo de Capacitancia
    charge = np.linalg.solve(mom, excv)
    charge_tot = np.sum(charge)
    cap_norm = charge_tot / unorm
    print(f"Capacitancia normalizada: {cap_norm:.4f}")

    # --- 5. GENERACIÓN DE GRÁFICAS ---
    # Densidad de carga: charge_celda / área_celda
    rho_s = charge / (a * b)
    # Redimensionar a matriz 2D para graficar
    rho_mat = rho_s.reshape((N, M), order='F')

    fig = plt.figure(figsize=(10, 7))
    ax = fig.add_subplot(111, projection='3d')
    
    surf = ax.plot_surface(xmat, ymat, rho_mat, cmap='viridis', edgecolor='black', lw=0.1)
    
    ax.set_title("Distribución de Densidad de Carga")
    ax.set_xlabel('x (m)')
    ax.set_ylabel('y (m)')
    ax.set_zlabel('Densidad de Carga Normalizada')
    fig.colorbar(surf, shrink=0.5, aspect=10)
    
    plt.show()

if __name__ == "__main__":
    solve_plate_vectorized()