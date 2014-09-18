COMPILEFLAGS = -O3
#LIBDIR = -L/usr/libs
LIBDIR = 
#INCLUDES = -I/usr/include
INCLUDES = 
LINKERFLAGS = -larmadillo

PROJECT = 2_program.x

obj = 2_program.o jacobi.o 
COMPILER = g++

default: $(PROJECT)

$(PROJECT): $(obj)
	$(COMPILER) $(LIBDIR) $(INCLUDES) -o $(PROJECT) $(obj) $(LINKERFLAGS)  
