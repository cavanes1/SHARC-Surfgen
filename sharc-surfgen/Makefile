#Macro definition
FC = ifort
FFLAGS =  -O3 -qopenmp -g -traceback -I$(MKLROOT)/include/intel64/lp64 -I$(MKLROOT)/include
LDFLAGS =  $(SURFGEN)/lib/libsurfgen.a \
					 $(MKLROOT)/lib/intel64/libmkl_blas95_lp64.a \
           $(MKLROOT)/lib/intel64/libmkl_lapack95_lp64.a \
           -Wl,--start-group \
           $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a \
           $(MKLROOT)/lib/intel64/libmkl_core.a \
           $(MKLROOT)/lib/intel64/libmkl_intel_thread.a \
           -Wl,--end-group -lpthread -lm
OBJSURFGEN= nn.o orient.o sub.o dms.o shepard.o pes.o soc.o
#end of Macro definition

sharc-surfgen.x:      $(OBJSURFGEN) sharc-surfgen.o
	$(FC) $(FFLAGS) $(OBJSURFGEN) sharc-surfgen.o -o sharc-surfgen.x $(LDFLAGS)

clean:
	rm -f *.o *.mod *.exe *.x a.out

%.o:	%.f90
	$(FC) $(FFLAGS) -c $< -o $@
