#arguments pour l'étape de compilation: on met tous les avertissements
CPPFLAGS := -Wall -Wextra

# arguments pour l'édition de liens, bibliothèques externes
LDFLAGS := -lm 


# nom du compilateur 
# CPP := g++ -I"C:\Users\remyl\eigen-3.4.0\Eigen"  
#Pas nécessaire si l'include path est ajouté dans "c-cpp_properties.json" 
CPP := g++ -I"C:/Users/remyl/eigen-3.4.0/Eigen" 

# Quelques remarques : -E fait uniquement la précompilation; -o <file> place l'output dans le fichier file; -c compile et assemble mais ne fait pas l'édition de liens


#cible par défaut
all: main


# compilation finale (avec édition de liens)
main: main.o
	$(CPP) -o main $^ $(LDFLAGS)
	./main


# compilation séparée
main.o: main.cpp
	$(CPP) -o $@ -c $< $(CPPFLAGS)


# Certains éditeurs créent aussi des .bak ou ~, mais apparemment pas ici 
clean: 
	rm -f *.o 


mrproper: clean 
	rm -f main

#aussi makedepend mais pas installé sur les ordis du magistère




