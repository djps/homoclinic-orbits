#
#   Makefile
#
PGM    = ext
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
run: first second third
#
first:	$(EXE)
	@echo " "
	@echo "Program $(PGM) is started"
	@echo " "
	@cp $(AUTO_DIR)/include/fcon.h fcon.h
	@echo "Starting conversion of $(PGM).dat : "
	@echo "(Required files : $(PGM).f, r.$(PGM), $(PGM).dat)"
	@f77 $(AUTO_DIR)/src/fcon.f $(PGM).f -o fcon $(AUTO_DIR)/lib/eispack.o
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
	@cp fort.7 p.$(PGM).1
	@cp fort.8 q.$(PGM).1
	@cp fort.9 d.$(PGM).1
	@cp q.$(PGM).1 t.$(PGM).1
	@echo " Saved as *.$(PGM)"
	@$(RM) fort.*
	@echo " Appended to *.$(PGM)"
#
second:	$(EXE)
	@echo " "
	@echo "$(PGM) : second run"
	@cp r.$(PGM).2 fort.2
	@cp q.$(PGM).1 fort.3
	@$(EXE)
	@cp fort.7 p.$(PGM).2
	@cp fort.8 q.$(PGM).2
	@cp fort.9 d.$(PGM).2
	@cp q.$(PGM).2 t.$(PGM).2
	@echo " Appended to *.$(PGM)"
	@$(RM) fort.*
#
third:	$(EXE)
	@echo " "
	@echo "$(PGM) : third run"
	@cp r.$(PGM).3 fort.2
	@cp q.$(PGM).2 fort.3
	@$(EXE)
	@cp fort.7 p.$(PGM).3
	@cp fort.8 q.$(PGM).3
	@cp fort.9 d.$(PGM).3
	@cp q.$(PGM).3 t.$(PGM).3
	@echo " Appended to *.$(PGM)"
	@$(RM) fort.*
#
fourth:	$(EXE)
	@echo " "
	@echo "$(PGM) : fourth run"
	@cp r.$(PGM).4 fort.2
	@cp t.$(PGM)   fort.3
	@$(EXE)
	@cat fort.7 >> p.$(PGM)
	@cat fort.8 >> q.$(PGM)
	@cat fort.9 >> d.$(PGM)
	@cp q.$(PGM) t.$(PGM).1
	@echo " Appended to *.$(PGM)"
	@$(RM) fort.*
#
fifth:	$(EXE)
	@echo " "
	@echo "$(PGM) : fifth run"
	@cp r.$(PGM).5 fort.2
	@cp t.$(PGM).1 fort.3
	@$(EXE)
	@cp fort.7 p.$(PGM).5
	@cp fort.8 q.$(PGM).5
	@cp fort.9 d.$(PGM).5
	@echo " Appended to *.$(PGM)"
	@$(RM) fort.*
#
sixth:	$(EXE)
	@echo " "
	@echo "$(PGM) : sixth run"
	@cp r.$(PGM).6 fort.2
	@cp t.$(PGM).1 fort.3
	@$(EXE)
	@cat fort.7 >> p.$(PGM).5
	@cat fort.8 >> q.$(PGM).5
	@cat fort.9 >> d.$(PGM).5
	@echo " Appended to *.$(PGM)"
	@$(RM) fort.*
#
seventh:$(EXE)
	@echo " "
	@echo "$(PGM) : seventh run"
	@cp r.$(PGM).7 fort.2
	@cp t.$(PGM).1 fort.3
	@$(EXE)
	@cp fort.7 p.$(PGM).7
	@cp fort.8 q.$(PGM).7
	@cp fort.9 d.$(PGM).7
	@echo " Appended to *.$(PGM)"
	@$(RM) fort.*
#
eighth:$(EXE)
	@echo " "
	@echo "$(PGM) : eighth run"
	@cp r.$(PGM).8 fort.2
	@cp t.$(PGM).1 fort.3
	@$(EXE)
	@cat fort.7 >> p.$(PGM).7
	@cat fort.8 >> q.$(PGM).7
	@cat fort.9 >> d.$(PGM).7
	@echo " Appended to *.$(PGM)"
	@$(RM) fort.*
#
ninth:$(EXE)
	@echo " "
	@echo "$(PGM) : ninth run"
	@cp r.$(PGM).9 fort.2
	@cp t.$(PGM).1 fort.3
	@$(EXE)
	@cp fort.7 p.$(PGM).9
	@cp fort.8 q.$(PGM).9
	@cp fort.9 d.$(PGM).9
	@echo " Appended to *.$(PGM)"
	@$(RM) fort.*
#
tenth:$(EXE)
	@echo " "
	@echo "$(PGM) : tenth run"
	@cp r.$(PGM).10 fort.2
	@cp t.$(PGM).1  fort.3
	@$(EXE)
	@cat fort.7 >> p.$(PGM).9
	@cat fort.8 >> q.$(PGM).9
	@cat fort.9 >> d.$(PGM).9
	@echo " Appended to *.$(PGM)"
	@$(RM) fort.*
#
clean:
	@echo "Cleaning $(PGM) ..."
	@cp $(SRC) $(EXE)
	@$(RM) *.exe *.o p.* q.* d.* fort.* *~ core
	@cp r.$(PGM).1 r.$(PGM)
	@echo "Cleaning ... done"
