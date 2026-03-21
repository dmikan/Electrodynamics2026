# Include the official MFEM configuration from the container
include /home/euler/mfem/config/config.mk

# Path to the MFEM directory (where mfem.hpp lives)
MFEM_DIR = /home/euler/mfem

# Search for all .cpp files in the current folder
SOURCES = $(wildcard *.cpp)

# Create a list of executables (removing the .cpp extension from each source)
PROGRAMS = $(SOURCES:.cpp=)

.SUFFIXES: .cpp .o
.PHONY: all clean

# By default, if you only type 'make', it tries to compile all .cpp files it finds
all: $(PROGRAMS)

# This rule says: "To create any executable '%', look for a '%' with .cpp extension"
%: %.cpp
	$(MFEM_CXX) $(MFEM_FLAGS) -I$(MFEM_DIR) $(MFEM_INC) $< -o $@ $(MFEM_LIBS)
	
# Clean up object files and executables
clean:
	rm -f *.o $(PROGRAMS)