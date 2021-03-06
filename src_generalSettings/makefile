#
# Makefile for 'createGeneralSettings'.
#
# Type 'make' or 'make newran' to create the binary.
# Type 'make clear' to delete all temporaries or 'make clean' to delete all temporaries and executables.
# Type 'make run' to execute the binary.
# Type 'make debug' to debug the binary using gdb(1).
#

# build target specs
CC = g++
CFLAGS = -O3 -std=c++11  
OUT_DIR = .
EXC_DIR = ..
LIBS =

# g++  maths.cpp newran02/newran.cpp newran02/extreal.cpp  newran02/myexcept.cpp main.cpp generateFiles.cpp -o test

# cible: dependance
#   commandes
# $@ => nom de la cible
# $< => nom de la 1ère dépendance
# $^ => liste des dépendances
# $? => liste des dépendances plus récentes que la cible
# $* => nom du fichier sans suffixe

# first target entry is the target invoked when typing 'make' ("default" or "all" are the most commonly used names by convention)
default: $(EXC_DIR)/createGeneralSettings

$(EXC_DIR)/createGeneralSettings: $(OUT_DIR)/main.cpp.o $(OUT_DIR)/maths.cpp.o $(OUT_DIR)/newran02/extreal.cpp.o $(OUT_DIR)/newran02/hist.cpp.o $(OUT_DIR)/newran02/myexcept.cpp.o $(OUT_DIR)/newran02/newran.cpp.o $(LIBS) $(OUT_DIR)/generateFiles.cpp.o
	@echo -n 'Linking createGeneralSettings... '
	@$(CC) $(CFLAGS) -o $@ $^
	@echo Done.

$(OUT_DIR)/maths.cpp.o: maths.cpp maths.h \
 newran02/newran.h newran02/include.h newran02/boolean.h \
 newran02/myexcept.h newran02/extreal.h
	@echo -n 'Compiling maths.cpp... '
	@$(CC) $(CFLAGS) -o $@ -c $<
	@echo Done.
	
$(OUT_DIR)/main.cpp.o: main.cpp newran02/newran.h \
 newran02/include.h newran02/boolean.h newran02/myexcept.h \
 newran02/extreal.h maths.h generateFiles.h
	@echo -n 'Compiling main.cpp... '
	@$(CC) $(CFLAGS) -o $@ -c $<
	@echo Done.
	
$(OUT_DIR)/generateFiles.cpp.o: generateFiles.cpp newran02/newran.h \
 newran02/include.h newran02/boolean.h newran02/myexcept.h \
 newran02/extreal.h maths.h generateFiles.h
	@echo -n 'Compiling generateFiles.cpp... '
	@$(CC) $(CFLAGS) -o $@ -c $<
	@echo Done.
	
$(OUT_DIR)/newran02/extreal.cpp.o: newran02/extreal.cpp newran02/include.h newran02/boolean.h newran02/extreal.h
	@echo -n 'Compiling extreal.cpp... '
	@$(CC) $(CFLAGS) -o $@ -c $<
	@echo Done.

$(OUT_DIR)/newran02/hist.cpp.o: newran02/hist.cpp newran02/include.h newran02/boolean.h newran02/extreal.h newran02/newran.h \
 newran02/myexcept.h newran02/tryrand.h
	@echo -n 'Compiling hist.cpp... '
	@$(CC) $(CFLAGS) -o $@ -c $<
	@echo Done.

$(OUT_DIR)/newran02/myexcept.cpp.o: newran02/myexcept.cpp newran02/include.h newran02/boolean.h newran02/myexcept.h
	@echo -n 'Compiling myexcept.cpp... '
	@$(CC) $(CFLAGS) -o $@ -c $<
	@echo Done.

$(OUT_DIR)/newran02/newran.cpp.o: newran02/newran.cpp newran02/include.h newran02/newran.h newran02/boolean.h \
 newran02/myexcept.h newran02/extreal.h
	@echo -n 'Compiling newran.cpp... '
	@$(CC) $(CFLAGS) -o $@ -c $< -w # add -w -> no warning (a lot because use deprecated conversion string -> char*)
	@echo Done.
	
run:
	./createGeneralSettings

debug:
	gdb ./createGeneralSettings
	
clear:
	@echo -n 'Removing all temporary binaries... '
	@rm -f $(OUT_DIR)/*.o  $(OUT_DIR)/newran02/*.o 
	@echo Done.

clean:
	@echo -n 'Removing executables and all temporary binaries... '
	@rm -f createGeneralSettings $(OUT_DIR)/*.o $(OUT_DIR)/newran02/*.o 
	@echo Done.

