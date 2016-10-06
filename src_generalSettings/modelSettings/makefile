#
# Makefile for 'createModelIN'.
#
# Type 'make' or 'make newran' to create the binary.
# Type 'make clean' or 'make clear' to delete all temporaries.
# Type 'make run' to execute the binary.
# Type 'make debug' to debug the binary using gdb(1).
#

# build target specs
CC = g++
CFLAGS = -O3 
OUT_DIR = .
LIBS =


#  g++ main.cpp maths.cpp newran02/newran.cpp newran02/extreal.cpp newran02/myexcept.cpp -o test


# createModelIN: main.cpp maths.h newran02/include.h newran02/boolean.h newran02/myexcept.h newran02/newran.h newran02/extreal.h
#	g++ -o createModelIN main.cpp maths.cpp newran02/newran.cpp newran02/extreal.cpp newran02/myexcept.cpp


# cible: dependance
#   commandes

# $@ => nom de la cible
# $< => nom de la 1ère dépendance
# $^ => liste des dépendances
# $? => liste des dépendances plus récentes que la cible
# $* => nom du fichier sans suffixe

# first target entry is the target invoked when typing 'make'
default: createModelIN

createModelIN: $(OUT_DIR)/main.cpp.o $(OUT_DIR)/maths.cpp.o $(OUT_DIR)/newran02/extreal.cpp.o $(OUT_DIR)/newran02/hist.cpp.o $(OUT_DIR)/newran02/myexcept.cpp.o $(OUT_DIR)/newran02/newran.cpp.o $(LIBS)
	@echo -n 'Linking createModelIN... '
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
 newran02/extreal.h maths.h settings.h
	@echo -n 'Compiling main.cpp... '
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
	@$(CC) $(CFLAGS) -o $@ -c $<
	@echo Done.
	
run:
	./createModelIN

debug:
	gdb ./createModelIN
	
clean:
	@echo -n 'Removing all temporary binaries... '
	@rm -f $(OUT_DIR)/*.o  $(OUT_DIR)/newran02/*.o 
	@echo Done.

clear:
	@echo -n 'Removing executables and all temporary binaries... '
	@rm -f createModelIN $(OUT_DIR)/*.o
	@echo Done.

