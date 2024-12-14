# arguments for compilation: we make the compilation verbose, and we ask for all warnings
CPPFLAGS := -Wall -Wextra

# arguments for linking 
LDFLAGS := -lm 


# Compiler and include path 
CPP := g++ -I"C:/Users/remyl/eigen-3.4.0/Eigen" 

# Default target 
all: main



# Compilation and linking 
main: main.o lattice_xy.o lattice_bh.o
	$(CPP) -o main $^ $(LDFLAGS)


# Separate compilation, so that not everything is recompiled when only one file changes
main.o: main.cpp src/lattice_xy.h src/lattice_bh.h
	$(CPP) -o $@ -c $< $(CPPFLAGS)


lattice_xy.o: src/lattice_xy.cpp src/lattice_xy.h
	$(CPP) -o $@ -c $< $(CPPFLAGS)


lattice_bh.o: src/lattice_bh.cpp src/lattice_bh.h
	$(CPP) -o $@ -c $< $(CPPFLAGS)



#----- Non-default target -----


# Target to run main.cpp
run: main 
	./main

# Target for tests and proof of concepts: it compiles and run directly the poc.cpp file 
poc: poc.o lattice_xy.o lattice_bh.o
	$(CPP) -o poc $^ $(LDFLAGS)
	./poc

poc.o: poc.cpp src/lattice_xy.h src/lattice_bh.h
	$(CPP) -o $@ -c $< $(CPPFLAGS)



# #Dependecy generation
# depend: 
# 	makedepend main.cpp

# # Include the dependencies		
# -include Makefile.dependencies

# # Rule to generate the dependency file
# Makefile.dependencies: main.cpp
# 	makedepend -f- main.cpp > Makefile.dependencies 2>/dev/null


# Clean target to remove object files
clean:
	rm -f *.o

# Mrproper target to remove object files and the executable
mrproper: clean
	rm -f *.exe


