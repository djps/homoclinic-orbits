#
#   Makefile
#
PGM    = mag
FC     = gfortran
FFLAGS = -O
SRC    = $(PGM).f
EXE    = $(PGM).exe
LIBS   = $(AUTO_DIR)/lib/*.o
RM     = rm -f
#
all: clean $(EXE) run
#
$(EXE): $(SRC:.f=.o)
	$(FC) $(FFLAGS) $(SRC:.f=.o) -o $@ $(LIBS)
#
run: first
#
first:	$(EXE)
	@echo " "
	@echo "Program $(PGM) is started"
	@echo " "
	@cp $(AUTO_DIR)/include/fcon.h fcon.h
	@echo "Starting conversion of $(PGM).dat : "
	@echo "(Required files : $(PGM).f, r.$(PGM), $(PGM).dat)"
	@f77 $(AUTO_DIR)/src/fcon.f $(PGM).f -o fcon $(AUTO_DIR)/lib/nag.o $(AUTO_DIR)/lib/eispack.o
	@cp r.$(PGM)   fort.2
	@cp $(PGM).dat fort.3
	@fcon
	@mv fort.8 q.dat
	@rm fcon* fort.2 fort.3
	@echo " Conversion done : converted file saved as q.dat"
	@echo " "
	@echo "$(PGM) : first run"
	@cp r.$(PGM).1 fort.2
	@cp q.dat fort.3
	@touch fort.3
	@cp $(PGM).dat fort.4
	@$(EXE)
	@cp fort.7 p.$(PGM)
	@cp fort.8 q.$(PGM)
	@cp fort.9 d.$(PGM)
	@cp fort.8 t.$(PGM)
	@echo " "
	@echo " Appended to *.$(PGM)"
#
second:	$(EXE)
	@echo " "
	@echo "$(PGM) : second run"
	@cp r.$(PGM).2 fort.2
	@cp t.$(PGM)   fort.3
	@$(EXE)	
	@cp fort.7 p.$(PGM).2
	@cp fort.8 q.$(PGM).2
	@cp fort.9 d.$(PGM).2
	@echo " Appended to *.$(PGM)"
	@$(RM) fort.*
	@echo " "
#
third:	$(EXE)
	@echo " "
	@echo "$(PGM) : third run "
	@cp r.$(PGM).3 fort.2
	@cp t.$(PGM)   fort.3
	@$(EXE)
	@cp fort.7 p.$(PGM).3
	@cp fort.8 q.$(PGM).3
	@cp fort.9 d.$(PGM).3
	@echo " Appended to *.$(PGM)"
	@echo " "
#
fourth:	$(EXE)
	@echo " "
	@echo "$(PGM) : forth run "
	@cp r.$(PGM).4 fort.2
	@cp q.$(PGM)   fort.3
	@$(EXE)
	@cat fort.7 >> p.$(PGM)
	@cat fort.8 >> q.$(PGM)
	@cat fort.9 >> d.$(PGM)
	@echo " Appended to *.$(PGM)"
	@echo " "
#
fifth:	$(EXE)
	@echo " "
	@echo "$(PGM) : fifth run : Limit Point Continuation with kappa, rho"
	@cp r.$(PGM).5 fort.2
	@cp q.$(PGM)   fort.3
	@$(EXE)
	@cp fort.7 p.$(PGM).2
	@cp fort.8 q.$(PGM).2
	@cp fort.9 d.$(PGM).2
	@echo " Appended to *.$(PGM)"
	@echo " "
#
clean:
	@echo "Cleaning $(PGM) ..."
	@cp $(SRC) $(EXE)
#	@$(RM) *.exe *.o p.* q.* d.* fort.* *~ core
	@$(RM) *.exe *.o p.* q.* d.* *~ core
	@cp r.$(PGM).1 r.$(PGM)
	@echo "Cleaning ... done"
