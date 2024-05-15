# Compilador y opciones
FC = gfortran
FFLAGS=-O2 -Wall -Wextra -std=legacy
# FFLAGS = -Wall -Wextra -g
# FFLAGS=-O2 -std=f2008
# FFLAGS=-O3 -Wall -Wextra

# Directorios
SRCDIR = .
MODDIR = modules
BINDIR = .

# Objetivos
SRCS := $(wildcard $(SRCDIR)/*.f90)
MODS := $(wildcard $(MODDIR)/*.f90)
OBJS := $(SRCS:.f90=.o) $(MODS:.f90=.o)
TARGET = $(BINDIR)/main

# Reglas de compilación
all: $(TARGET)

$(TARGET): $(OBJS)
	$(FC) $(FFLAGS) -o $@ $(OBJS)

$(SRCDIR)/%.o: $(SRCDIR)/%.f90
	$(FC) $(FFLAGS) -c -o $@ $<

$(MODDIR)/%.o: $(MODDIR)/%.f90
	$(FC) $(FFLAGS) -c -o $@ $<

# Limpieza
clean:
	rm -f $(TARGET) $(OBJS)
