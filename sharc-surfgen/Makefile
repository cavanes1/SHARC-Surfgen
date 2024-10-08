#Macro definition
FC = ifort
FFLAGS:=-O3 -qopenmp -I$(MKLROOT)/include/intel64/lp64 -I$(MKLROOT)/include
LDFLAGS =  $(SURFGEN)/lib/libsurfgen.a \
					 $(MKLROOT)/lib/intel64/libmkl_blas95_lp64.a \
           $(MKLROOT)/lib/intel64/libmkl_lapack95_lp64.a \
           -Wl,--start-group \
           $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a \
           $(MKLROOT)/lib/intel64/libmkl_core.a \
           $(MKLROOT)/lib/intel64/libmkl_intel_thread.a \
           -Wl,--end-group -lpthread -lm
LIBSURF =  hddata.o combinatorial.o progdata.o libutil.o libsym.o libinternal.o \
           io.o potlib.o getver.o
OBJDIP = nn.o dms.o shepard.o sub.o orient.o dip.o
OBJSOC = nn.o soc.o sub.o orient.o qm.o
#end of Macro definition

all:	lib surfgen dip qm
	
surfgen:
	$(FC) $(FFLAGS) -c surfgen.f90
	$(FC) $(FFLAGS) -o surfgen.x surfgen.o libsurfgen.a $(LDFLAGS)

dip:	$(OBJDIP)
	$(FC) $(FFLAGS) $(OBJDIP) -o dip.x $(LDFLAGS)

qm:	$(OBJSOC)
	$(FC) $(FFLAGS) $(OBJSOC) -o qm.x $(LDFLAGS)
	rm -f *.o *.mod

lib:	$(LIBSURF)
	$(AR) -r -v libsurfgen.a $(LIBSURF)
	@echo " Cleaning up..."
	rm -f $(LIBSURF)
	rm -f ./*.mod

./%.o : ./lib/%.f90
	@echo ' Building file: $<'
	$(FC) $(FCOPT) -c -o $@ $<
	@echo ' Finished building $<'
	@echo ''

./%.o : ./lib/%.F90
	@echo ' Building file: $<'
	$(FC) $(FCOPT) -c -o $@ $<
	@echo ' Finished building $<'
	@echo ''

./%.o : ./%.f90
	@echo ' Building file: $<'
	$(FC) $(FCOPT) -c -o $@ $<
	@echo ' Finished building $<'
	@echo ''

clean:
	rm -f *.o *.mod *.exe *.x a.out

%.o:	%.f90
	$(FC) $(FFLAGS) -c $< -o $@
