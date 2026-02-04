all : clean shooting

FFLAGS = -ffixed-line-length-132

LIB_LIST = -llapack -lblas
#
shooting.o : shooting.f
	gfortran  $(FFLAGS) -c shooting.f 
#
monodromy.o : monodromy.f
	gfortran  $(FFLAGS) -c monodromy.f
#
newton.o : newton.f
	gfortran  $(FFLAGS) -c newton.f 
#
postprocess.o : postprocess.f
	gfortran  $(FFLAGS) -c postprocess.f
#
dop853.o : dop853.f
	gfortran  $(FFLAGS) -c dop853.f 
#
shooting : shooting.o monodromy.o newton.o postprocess.o dop853.o
	gfortran  -o run shooting.o monodromy.o newton.o postprocess.o dop853.o $(LIB_LIST)
#
clean: 
	@echo "Cleaning program ..."
	rm -rf *.o core *.exe *~ run
	@echo "Cleaning ... done"