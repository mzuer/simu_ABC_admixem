#
# Makefile for 'runABC'.
#
# Type 'make' or 'make newran' to create the binary.
# Type 'make clear' to delete all temporaries or 'make clean' to delete all temporaries and executables.
# Type 'make run' to execute the binary.
# Type 'make debug' to debug the binary using gdb(1).
#

# g++  -std=c++11 `pkg-config --cflags --libs gsl` -lboost_regex

# build target specs
CC = g++
CFLAGS = -O3 -std=c++11 -fopenmp
OUT_DIR = src
SRC_DIR = src
EXEC_DIR = .
LIBS =
#GSLFLAGS = -I/usr/local/include  -L/usr/local/lib -lgsl -lgslcblas #  pkg-config --cflags --libs gsl
GSLFLAGS = `pkg-config --cflags --libs gsl`
BOOSTFLAGS = -lboost_regex

# cible: dependance
#   commandes
# $@ => nom de la cible
# $< => nom de la 1ère dépendance
# $^ => liste des dépendances
# $? => liste des dépendances plus récentes que la cible
# $* => nom du fichier sans suffixe

# first target entry is the target invoked when typing 'make' ("default" or "all" are the most commonly used names by convention)
default: runABC

$(EXEC_DIR)/runABC: $(OUT_DIR)/main.cpp.o $(OUT_DIR)/model.cpp.o $(OUT_DIR)/simu.cpp.o $(OUT_DIR)/my_math.cpp.o $(OUT_DIR)/summarystat.cpp.o $(OUT_DIR)/retrieve_input.cpp.o
	@echo -n 'Linking runABC... '
	@$(CC) $(CFLAGS) -o $@ $^ $(GSLFLAGS) $(BOOSTFLAGS)
	@echo Done.

$(OUT_DIR)/main.cpp.o: $(SRC_DIR)/main.cpp $(SRC_DIR)/simu.h $(SRC_DIR)/model.h $(SRC_DIR)/retrieve_input.h
	@echo -n 'Compiling main.cpp... '
	@$(CC) $(CFLAGS) -o $@ -c $< $(GSLFLAGS) $(BOOSTFLAGS)
	@echo Done.
	
$(OUT_DIR)/model.cpp.o: $(SRC_DIR)/model.cpp $(SRC_DIR)/settings.h $(SRC_DIR)/model.h $(SRC_DIR)/my_math.h
	@echo -n 'Compiling model.cpp... '
	@$(CC) $(CFLAGS) -o $@ -c $< $(GSLFLAGS) $(BOOSTFLAGS)
	@echo Done.
	
$(OUT_DIR)/simu.cpp.o: $(SRC_DIR)/simu.cpp $(SRC_DIR)/settings.h $(SRC_DIR)/simu.h $(SRC_DIR)/my_math.h $(SRC_DIR)/summarystat.h
	@echo -n 'Compiling simu.cpp... '
	@$(CC) $(CFLAGS) -o $@ -c $< $(GSLFLAGS) $(BOOSTFLAGS)
	@echo Done.
	
$(OUT_DIR)/summarystat.cpp.o: $(SRC_DIR)/summarystat.cpp $(SRC_DIR)/summarystat.h $(SRC_DIR)/settings.h
	@echo -n 'Compiling summarystat.cpp... '
	@$(CC) $(CFLAGS) -o $@ -c $< $(GSLFLAGS) $(BOOSTFLAGS)
	@echo Done.
	
$(OUT_DIR)/my_math.cpp.o: $(SRC_DIR)/my_math.cpp $(SRC_DIR)/my_math.h
	@echo -n 'Compiling my_math.cpp... '
	@$(CC) $(CFLAGS) -o $@ -c $< $(GSLFLAGS) $(BOOSTFLAGS)       #$(GSLFLAGS) $(BOOSTFLAGS) # BOOST for the regex
	@echo Done.
		
	
$(OUT_DIR)/retrieve_input.cpp.o: $(SRC_DIR)/retrieve_input.cpp $(SRC_DIR)/retrieve_input.h
	@echo -n 'Compiling retrieve_input.cpp... '
	@$(CC) $(CFLAGS) -o $@ -c $< #$(GSLFLAGS) $(BOOSTFLAGS)       #$(GSLFLAGS) $(BOOSTFLAGS) # BOOST for the regex
	@echo Done.
				
		
		
run:
	./runABC

debug:
	gdb ./runABC
	
clear:
	@echo -n 'Removing all temporary binaries... '
	@rm -f $(OUT_DIR)/*.o  
	@echo Done.

clean:
	@echo -n 'Removing executables and all temporary binaries... '
	@rm -f $(EXEC_DIR)/runABC $(OUT_DIR)/*.o 
	@echo Done.

