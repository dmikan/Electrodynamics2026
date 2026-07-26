# 1. Configuración de MFEM
include /home/euler/mfem/config/config.mk
MFEM_DIR = /home/euler/mfem

# 2. Directorios de tu proyecto
INCDIR = include
SRCDIR = src
APPDIR = apps
BUILDDIR = build

# 3. Compilador y banderas
# Usamos las banderas de MFEM y añadimos nuestra carpeta 'include'
CXX = $(MFEM_CXX)
CXXFLAGS = $(MFEM_FLAGS) -I$(MFEM_DIR) -I$(INCDIR) $(MFEM_INC)
LIBS = $(MFEM_LIBS)

# 4. Archivos
# Buscamos todos los .cpp en src/ (tus clases)
SOURCES = $(wildcard $(SRCDIR)/*.cpp)
# Convertimos src/Clase.cpp en build/Clase.o
OBJECTS = $(SOURCES:$(SRCDIR)/%.cpp=$(BUILDDIR)/%.o)

# Nombre del ejecutable final
TARGET = poisson_delta

# --- REGLAS DE COMPILACIÓN ---

all: $(BUILDDIR) $(TARGET)

# Crear la carpeta build si no existe
$(BUILDDIR):
	mkdir -p $(BUILDDIR)

# Regla para el ejecutable final (Enlace)
$(TARGET): $(APPDIR)/main.cpp $(OBJECTS)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LIBS)

# Regla para compilar cada clase (.cpp -> .o)
$(BUILDDIR)/%.o: $(SRCDIR)/%.cpp $(INCDIR)/%.hpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -rf $(BUILDDIR) $(TARGET)
	rm -f solution.gf mesh.mesh

.PHONY: all clean