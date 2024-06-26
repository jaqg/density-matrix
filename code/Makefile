# Variables
FC = gfortran
FCFLAGS = -g -Wall -O2
TARGET = main
MODULES_DIR = modules
SOURCES = $(wildcard $(MODULES_DIR)/*.f90) main.f90
OBJECTS = $(SOURCES:.f90=.o)
LDFLAGS = -llapack -lblas

# Default target
all: $(TARGET)

# Rule to build the target
$(TARGET): $(OBJECTS)
	$(FC) $(FCFLAGS) -o $@ $(OBJECTS) $(LDFLAGS)

# Rule to build object files from source files
%.o: %.f90
	$(FC) $(FCFLAGS) -c $< -o $@

test1:$(TARGET)
	@cp $(TARGET) test1
	@cd test1
	@./$(TARGET)

# Clean target
clean:
	rm -f $(OBJECTS) $(TARGET)

# Phony targets
.PHONY: all clean test1
