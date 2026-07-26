import numpy as np
import matplotlib.pyplot as plt

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
                         + (v**2 + zc**2)*np.arcsinh(u/np.sqrt(v**2 + zc**2))))
    prim4 = lambda u, v: (1/3*(u**2 + v**2 + zc**2)**(3/2))
    
    # Cálculo de las contribuciones en los diferentes límites (Yuu, Yud, Ydu, Ydd)
    Yuu = ((xc - a)*(yc - b)*definite(prim1, xc-a, xc, yc-b, yc)
           - (yc - b)*definite(prim2, xc-a, xc, yc-b, yc)
           - (xc - a)*definite(prim3, xc-a, xc, yc-b, yc)
           + definite(prim4, xc-a, xc, yc-b, yc))
    Yud = (-(xc - a)*(yc + b)*definite(prim1, xc-a, xc, yc, yc+b)
           + (yc + b)*definite(prim2, xc-a, xc, yc, yc+b)
           + (xc - a)*definite(prim3, xc-a, xc, yc, yc+b)
           - definite(prim4, xc-a, xc, yc, yc+b))
    Ydu = (-(xc + a)*(yc - b)*definite(prim1, xc, xc+a, yc-b, yc)
           + (yc - b)*definite(prim2, xc, xc+a, yc-b, yc)
           + (xc + a)*definite(prim3, xc, xc+a, yc-b, yc)
           - definite(prim4, xc, xc+a, yc-b, yc))
    Ydd = ((xc + a)*(yc + b)*definite(prim1, xc, xc+a, yc, yc+b)
           - (yc + b)*definite(prim2, xc, xc+a, yc, yc+b)
           - (xc + a)*definite(prim3, xc, xc+a, yc, yc+b)
           + definite(prim4, xc, xc+a, yc, yc+b))
    # El resultado final se normaliza por el producto de las áreas de ambas celdas
    return (Yuu + Yud + Ydu + Ydd) / (a * b)**2


# --- 2. SOLUCION PLACA ESTÁTICA (Resultados del articulo)---
def solve_plate_vectorized():
    # Configuración de la placa
    aside, bside = 1.0, 1.0
    M, N = 20, 20  # Divisiones por lado
    a, b = aside / M, bside / N
    tresh = 10     # Umbral de aproximación
    flag_g = 0     # 0 para Point-Matching, 1 para Galerkin

    # Mallado de centros de celdas
    x = np.linspace(a/2, aside - a/2, M)
    y = np.linspace(b/2, bside - b/2, N)
    xmat, ymat = np.meshgrid(x, y)
    
    # Aplanado de vectores
    xvec = xmat.flatten(order='F')
    yvec = ymat.flatten(order='F')

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
        metodo = "Galerkin"
    else:
        mom = np.where(rh > threshold_val, 1/R, DiscrGF(XC, YC, ZC, a, b))
        metodo = "Point-Matching"

    # Vector de excitacion
    
    # ESCENARIO 1: Potencial constante unitario
    #unorm = 1.0
    #excv = unorm * np.ones((ndim, 1))
    
    # ESCENARIO 2: Carga puntual unitaria
    xq, yq, zq = 0.5, 0.5, 0.5  # Posición de la carga puntual sobre el centro
    excv = 1.0 / np.sqrt((xq - xvec)**2 + (yq - yvec)**2 + zq**2)
    excv = excv.reshape(-1, 1)

    # Resolución y cálculo de Capacitancia
    charge = np.linalg.solve(mom, excv)
    charge_tot = np.sum(charge)
    #cap_norm = charge_tot / unorm
    #print(f"Capacitancia normalizada: {cap_norm:.4f}")
    print(f"Carga total: {charge_tot:.4f}")

    # -- GENERACIÓN DE GRÁFICAS ---
    # Densidad de carga: charge_celda / área_celda
    rho_s = charge / (a * b)
    # Redimensionar a matriz 2D para graficar
    rho_mat = rho_s.reshape((N, M), order='F')

    fig = plt.figure(figsize=(10, 7))
    ax = fig.add_subplot(111, projection='3d')
    
    surf = ax.plot_surface(xmat, ymat, rho_mat, cmap='viridis', edgecolor='black', 
                           lw=0.1)
    
    ax.set_title(f"Distribución de Densidad de Carga\n({metodo})")
    ax.set_xlabel('x (m)')
    ax.set_ylabel('y (m)')
    ax.set_zlabel('Densidad de Carga Normalizada')
    fig.colorbar(surf, shrink=0.5, aspect=10)
    
    plt.show()
    
# --- 3. CILINDRO CERRADO ---
def potential_cylinder(xvec, yvec, zvec, charge_int, charge_ext,
                       pos_q1, pos_q2, R):

    # Malla de observación (plano z = 0.5)
    xo = np.linspace(-1.5, 1.5, 100)
    yo = np.linspace(-1.5, 1.5, 100)

    Xo, Yo = np.meshgrid(xo, yo)
    Zo = 0.5

    fig, axes = plt.subplots(1, 2, figsize=(16, 7))

    casos = [
        (charge_int, pos_q1, "Carga interna (0,0,0.5)"),
        (charge_ext, pos_q2, "Carga externa (0.8,0,0.5)")
    ]

    for i, (charge, pos_q, titulo) in enumerate(casos):

        xq, yq, zq = pos_q

        # Potencial de la carga puntual
        dist_q = np.sqrt((Xo - xq)**2 + (Yo - yq)**2 + (Zo - zq)**2)

        V_q = 1.0 / (dist_q + 1e-6)

        # Potencial debido a las cargas inducidas
        dist_cells = np.sqrt((Xo[:, :, np.newaxis] - xvec)**2 + 
                             (Yo[:, :, np.newaxis] - yvec)**2 + (Zo - zvec)**2)

        V_ind = np.sum(charge / (dist_cells + 1e-6), axis=2)

        # Potencial total
        Potencial = V_q + V_ind

        # Gráfica
        ax = axes[i]

        niveles = np.linspace(-0.5, 0.5, 21)

        cp = ax.contourf(Xo, Yo, Potencial, levels=niveles, cmap='RdBu_r', 
                         extend='both')

        ax.contour(Xo, Yo, Potencial, levels=niveles, colors='black', 
                   linewidths=0.5, alpha=0.5)

        # Cilindro
        circle = plt.Circle((0, 0), R, fill=False, color='black', linestyle='--',
                            linewidth=2, label='Cilindro')

        ax.add_artist(circle)

        # Carga fuente
        ax.scatter(xq, yq, color='yellow', edgecolor='black', s=100, label='Fuente (q=1)', 
                   zorder=5)

        ax.set_title(titulo)
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_aspect('equal')

        ax.legend()

        fig.colorbar(cp, ax=ax, label=r'Potencial Normalizado ($4\pi\varepsilon_0V$)')

    plt.tight_layout()
    plt.show()

def solve_cylinder():
    # --- PARÁMETROS GEOMÉTRICOS ---
    R = 0.5   # Radio del cilindro
    H = 1.0   # Altura del cilindro
    M = 32    # Divisiones angulares (horizontales)
    N = 24    # Divisiones verticales
    
    a_cell = (2 * np.pi * R) / M  # Ancho de la celda (arco)
    b_cell = H / N                # Alto de la celda
    
    flag_g = 0 # 0 para Point-Matching, 1 para Galerkin
      
    # --- MALLADO ---
    theta = np.linspace(0, 2*np.pi, M, endpoint=False)
    z = np.linspace(b_cell/2, H - b_cell/2, N)
    theta_mat, z_mat = np.meshgrid(theta, z)
    
    # Coordenadas cartesianas de los centros
    xvec = (R * np.cos(theta_mat)).flatten()
    yvec = (R * np.sin(theta_mat)).flatten()
    zvec = z_mat.flatten()
    
    # --- LLENADO DE LA MATRIZ MoM ---
    dx = xvec[:, np.newaxis] - xvec[np.newaxis, :]
    dy = yvec[:, np.newaxis] - yvec[np.newaxis, :]
    dz = zvec[:, np.newaxis] - zvec[np.newaxis, :]
    R_dist = np.sqrt(dx**2 + dy**2 + dz**2)
    
    # Llenado con 1/r (GF aproximada) para interacciones lejanas
    mom = 1.0 / (R_dist + np.finfo(float).eps)
    
    if flag_g == 1:
        self_val = AvDiscrGF(0, 0, 0, a_cell, b_cell) # Galerkin
    else:
        self_val = DiscrGF(0, 0, 0, a_cell, b_cell)   # Point-Matching
        
    np.fill_diagonal(mom, self_val)
    
    # --- EXCITACIÓN: CARGA PUNTUAL INTERNA ---
    # Escenario 1: Carga puntual dentro del cilindro

    xq1, yq1, zq1 = 0.0, 0.0, 0.5

    dist_q1 = np.sqrt((xvec - xq1)**2 + (yvec - yq1)**2 + (zvec - zq1)**2)
    excv_int= 1.0 / dist_q1  # Vector de excitación (potencial inducido por q)
    
    # Escenario 2: Carga puntual fuera del cilindro
    xq2, yq2, zq2 = 0.8, -0.5, 0.5

    dist_q2 = np.sqrt((xvec - xq2)**2 + (yvec - yq2)**2 + (zvec - zq2)**2)
    excv_ext = 1.0 / dist_q2  # Vector de excitación (potencial inducido por q)
    
    # --- RESOLUCIÓN DEL SISTEMA ---
    charge_int = np.linalg.solve(mom, -excv_int)
    charge_ext = np.linalg.solve(mom, -excv_ext)
        
    # --- 3. VISUALIZACIÓN ---
    # Creamos las mallas de coordenadas para graficar
    X_plot = R * np.cos(theta_mat)
    Y_plot = R * np.sin(theta_mat)
    Z_plot = z_mat
    
    # AJUSTE GRÁFICO: Duplicamos la columna inicial al final para cerrar el cilindro
    X_plot = np.hstack((X_plot, X_plot[:, :1]))
    Y_plot = np.hstack((Y_plot, Y_plot[:, :1]))
    Z_plot = np.hstack((Z_plot, Z_plot[:, :1]))
    
    fig, axes = plt.subplots(1, 2, figsize=(16, 7), subplot_kw={'projection': '3d'})
    metodo = "Galerkin" if flag_g == 1 else "Point-Matching"
    
    # --- 4. GRÁFICA 3D ---
    for i, (charge, pos_q, title, cmap) in enumerate([
        (charge_int, (xq1,yq1,zq1), "Carga Interna", 'viridis_r'),
        (charge_ext, (xq2,yq2,zq2), "Carga Externa", 'viridis_r')
    ]):
        rho_mat = (charge / (a_cell * b_cell)).reshape((N, M))
        rho_v = np.hstack((rho_mat, rho_mat[:, :1]))
        
        ax = axes[i]
        norm = plt.Normalize(rho_v.min(), rho_v.max())
        surf = ax.plot_surface(X_plot, Y_plot, Z_plot, 
                               facecolors=plt.cm.get_cmap(cmap)(norm(rho_v)), 
                               shade=False, edgecolor='black', lw=0.1, alpha=0.8)
        ax.scatter(*pos_q, color='black', s=100, label=f'Fuente (q=1)')
        ax.set_title(f"{title}\n({metodo})")
        fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax, 
                     shrink=0.4, label='Densidad $\\rho_s$')
        print(f"{title}: Carga total inducida = {np.sum(charge):.4f}")
    
    plt.legend()
    plt.tight_layout()
    plt.show()
    
    potential_cylinder(xvec, yvec, zvec, charge_int, charge_ext, (xq1,yq1,zq1), (xq2,yq2,zq2), R)


#if __name__ == "__main__":
#    solve_plate_vectorized()
if __name__ == "__main__":
    solve_cylinder()