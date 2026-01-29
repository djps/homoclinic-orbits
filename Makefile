all : clean shooting

FFLAGS = -ffixed-line-length-132
#
shooting.o : shooting.f
	f77 $(FFLAGS) -c shooting.f 
#
monodromy.o : monodromy.f
	f77 $(FFLAGS) -c monodromy.f
#
newton.o : newton.f
	f77 $(FFLAGS) -c newton.f 
#
postprocess.o : postprocess.f
	f77 $(FFLAGS) -c postprocess.f
#
dop853.o : dop853.f
	f77 $(FFLAGS) -c dop853.f 
#
shooting : shooting.o monodromy.o newton.o postprocess.o dop853.o
	f77 -o run shooting.o monodromy.o newton.o postprocess.o dop853.o
#
clean: 
	@echo "Cleaning program ..."
	rm -rf *.o core *.exe *~ run
	@echo "Cleaning ... done"