# Variables
FC = gfortran
FCFLAGS = -g -Wall -O2
TARGET = main
MODULES_DIR = modules
SOURCES = $(wildcard $(MODULES_DIR)/*.f90) main.f90
OBJECTS = $(SOURCES:.f90=.o)

# Default target
all: $(TARGET)

# Rule to build the target
$(TARGET): $(OBJECTS)
	$(FC) $(FCFLAGS) -o $@ $(OBJECTS)

# Rule to build object files from source files
%.o: %.f90
	$(FC) $(FCFLAGS) -c $< -o $@

# Clean target
clean:
	rm -f $(OBJECTS) $(TARGET)

# Phony targets
.PHONY: all clean
