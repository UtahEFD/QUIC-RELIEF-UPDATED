PROG =	quicurb_LIN64.exe
#PROG =	quicurb_MACI.exe

SRCS =	bisect.f90 bridgewake.f90 building_connect.f90 \
        building_damage.f90 building_parameterizations.f90 canopy_slope_match.f90 courtyard.f90 cylinderwake.f90 datamodule.f90 \
        defbuild.f90 diffusion.f90 divergence.f90 euler.f90 init.f90 main.f90 outfile.f90 \
        parking_garage.f90 pentagon.f90 plantinit.f90 poisson.f90 polywake.f90\
        read_hotmac_met.f90 read_ittmm5_met.f90 read_quic_met.f90 \
        reliefwake.f90 rectanglewake.f90 regress.f90 reliefrooftop.f90 rooftop.f90 sensorinit.f90 sidewall.f90 sor3d.f90 \
        sort.f90 street_intersect.f90 streetcanyon.f90 surface_coords.f90 turbulence_model.f90 upwind.f90 \
        utmll.f90 wallbc.f90 zd_bisect.f90 \
            
OBJS =	bisect.o bridgewake.o building_connect.o \
        building_damage.o building_parameterizations.o canopy_slope_match.o courtyard.o cylinderwake.o datamodule.o \
        defbuild.o diffusion.o divergence.o euler.o init.o main.o outfile.o \
        parking_garage.o pentagon.o plantinit.o poisson.o polywake.o\
        read_hotmac_met.o read_ittmm5_met.o read_quic_met.o \
        reliefwake.o rectanglewake.o regress.o reliefrooftop.o rooftop.o sensorinit.o sidewall.o sor3d.o \
        sort.o street_intersect.o streetcanyon.o surface_coords.o turbulence_model.o upwind.o \
        utmll.o wallbc.o zd_bisect.o \


LIBS =	

export MACOSX_DEPLOYMENT_TARGET=10.10

CC = cc
CFLAGS = -O
# compile with absoft

#FC = f95
#FFLAGS = -lib -f -s -N15 -N27 -O2
#F95 = f95
#F95FLAGS = -N27 -O3 -OPT:Olimit=0 -m64 #-mcmodel=medium
#LDFLAGS = -m64

# Compile with PGI

#FC = pgf95
#FFLAGS = -fastsse
#F95 = pgf95
#F95FLAGS = -fastsse

# compile with gfortran

#FC = gfortran
#FFLAGS = -O2 -fdefault-double-8 -fdefault-real-8 -fdefault-integer-8 -ffree-form
F95 = gfortran
F95FLAGS = -O3 -m64 -ffree-form -fconvert=little-endian -fopenmp
LDFLAGS = -m64 -fopenmp

# compile with intel fortran

#FC = ifort
#FFLAGS = -O2 -r8 -i8 -autodouble -xW -free -convert little_endian -i-static
#F95 = ifort
#F95FLAGS = -O3 -free -m64 -convert little_endian -openmp -static-intel
#LDFLAGS = -m64 -openmp -static-intel

all: $(PROG)

$(PROG): $(OBJS)
	$(F95) $(LDFLAGS) -o $@ $(OBJS) $(LIBS)

#DATE := /bin/date
#TODAY := `$(DATE) +%d-%h-%y`
#TARFILE := quic-urb.$(TODAY).tar
#COMPRESS := gzip

#tar:
#	tar -cf $(TARFILE) Makefile $(SRCS)
#	$(COMPRESS) -v $(TARFILE)
#clean:
#	rm -f $(PROG) $(OBJS) *.mod

.SUFFIXES: $(SUFFIXES) .f90

.f90.o:
	$(F95) $(F95FLAGS) -c $<

building_connect.o: datamodule.o
building_parameterizations.o: datamodule.o
defbuild.o: datamodule.o
divergence.o: datamodule.o
euler.o: datamodule.o
init.o: datamodule.o
main.o: datamodule.o
outfile.o: datamodule.o
pentagon.o: datamodule.o
plantinit.o: datamodule.o
polywake.o: datamodule.o
regress.o: datamodule.o
sensorinit.o: datamodule.o
sor3d.o: datamodule.o
sort.o: datamodule.o
street_intersect.o: datamodule.o
upwind.o: datamodule.o
wallbc.o: datamodule.o
poisson.o: datamodule.o
read_quic_met.o: datamodule.o
read_ittmm5_met.o: datamodule.o
read_hotmac_met.o: datamodule.o
cylinderwake.o: datamodule.o
reliefwake.o: datamodule.o
reliefrooftop.o: datamodule.o
rectanglewake.o: datamodule.o
wake.o: datamodule.o
streetcanyon.o: datamodule.o
relief.o: datamodule.o
rooftoop.o: datamodule.o
parking_garage.o: datamodule.o
courtyard.o: datamodule.o
bridgewake.o: datamodule.o
building_damage.o: datamodule.o
surface_coords.o: datamodule.o
sidewall.o: datamodule.o
