import numpy as np

def definite(prim, x1, x2, y1, y2):
    """
    Calcula la integral doble de una función f(x,y) sobre un dominio
    rectangular [x1, x2] x [y1, y2] evaluando la primitiva en las 4 esquinas.
    """
    return prim(x1, y1) + prim(x2, y2) - prim(x1, y2) - prim(x2, y1)

def DiscrGF(xc, yc, zc, a, b):
    """
    Calcula la Función de Green Discreta para Point-Matching (Gamma).
    """
    def prim1(u, v, h):
        return ((u * np.arcsinh(v / np.sqrt(u**2 + h**2))) +
                (v * np.arcsinh(u / np.sqrt(v**2 + h**2))) -
                (h * np.arctan(u * v / (h * np.sqrt(u**2 + v**2 + h**2)))))

    zc_safe = zc + np.finfo(float).eps
    result = (prim1(xc + a/2, yc + b/2, zc_safe) +
              prim1(xc - a/2, yc - b/2, zc_safe) -
              prim1(xc + a/2, yc - b/2, zc_safe) -
              prim1(xc - a/2, yc + b/2, zc_safe))
    return result / (a * b)

def AvDiscrGF(xc, yc, zc, a, b):
    """
    Calcula la interacción promedio entre dos celdas rectangulares (Galerkin, Psi).
    Reduce la integral 4D a 2D utilizando 4 primitivas analíticas.
    """
    zc_safe = zc + np.finfo(float).eps

    prim1 = lambda u, v: ((u * np.arcsinh(v / np.sqrt(u**2 + zc_safe**2))) +
                          (v * np.arcsinh(u / np.sqrt(v**2 + zc_safe**2))) -
                          (zc_safe * np.arctan(u * v / (zc_safe * np.sqrt(u**2 + v**2 + zc_safe**2)))))

    prim2 = lambda u, v: (0.5 * (v * np.sqrt(u**2 + v**2 + zc_safe**2) +
                          (u**2 + zc_safe**2) * np.arcsinh(v / np.sqrt(u**2 + zc_safe**2))))

    prim3 = lambda u, v: (0.5 * (u * np.sqrt(u**2 + v**2 + zc_safe**2) +
                          (v**2 + zc_safe**2) * np.arcsinh(u / np.sqrt(v**2 + zc_safe**2))))

    prim4 = lambda u, v: ((1.0 / 3.0) * (u**2 + v**2 + zc_safe**2)**(1.5))

    Yuu = ((xc - a)*(yc - b)*definite(prim1, xc-a, xc, yc-b, yc) -
           (yc - b)*definite(prim2, xc-a, xc, yc-b, yc) -
           (xc - a)*definite(prim3, xc-a, xc, yc-b, yc) +
           definite(prim4, xc-a, xc, yc-b, yc))

    Yud = (-(xc - a)*(yc + b)*definite(prim1, xc-a, xc, yc, yc+b) +
           (yc + b)*definite(prim2, xc-a, xc, yc, yc+b) +
           (xc - a)*definite(prim3, xc-a, xc, yc, yc+b) -
           definite(prim4, xc-a, xc, yc, yc+b))

    Ydu = (-(xc + a)*(yc - b)*definite(prim1, xc, xc+a, yc-b, yc) +
           (yc - b)*definite(prim2, xc, xc+a, yc-b, yc) +
           (xc + a)*definite(prim3, xc, xc+a, yc-b, yc) -
           definite(prim4, xc, xc+a, yc-b, yc))

    Ydd = ((xc + a)*(yc + b)*definite(prim1, xc, xc+a, yc, yc+b) -
           (yc + b)*definite(prim2, xc, xc+a, yc, yc+b) -
           (xc + a)*definite(prim3, xc, xc+a, yc, yc+b) +
           definite(prim4, xc, xc+a, yc, yc+b))

    return (Yuu + Yud + Ydu + Ydd) / (a * b)**2
