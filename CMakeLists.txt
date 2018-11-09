CMAKE_MINIMUM_REQUIRED (VERSION 2.8) 
PROJECT (QUICUrb)

ENABLE_LANGUAGE (Fortran)

if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
   # Apple/OS X based systems
   SET(EXECARCH "_MACI.exe")
elseif(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
   # Linux based system
   set(EXECARCH "_LIN64.exe")
elseif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
   set(EXECARCH "_WIN64.exe")
else(UNKNOWN)
   set(EXECARCH "_UNKNOWNARCH.exe")
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

SET (APPNAME "quicurb${EXECARCH}")
MESSAGE(STATUS "Setting app name to ${APPNAME}")

ADD_EXECUTABLE( ${APPNAME}
		bisect.f90 
		bridgewake.f90 
		building_connect.f90 
		building_damage.f90
		building_parameterizations.f90
		courtyard.f90
		cylinderwake.f90
		datamodule.f90
		defbuild.f90 
		diffusion.f90 
		divergence.f90 
		euler.f90 
		init.f90 
		main.f90 
		outfile.f90
		parking_garage.f90 
		pentagon.f90 
		plantinit.f90 
		poisson.f90 
		polywake.f90
        read_hotmac_met.f90 
		read_ittmm5_met.f90 
		read_quic_met.f90
		rectanglewake.f90 
		regress.f90 
		reliefwake.f90 
		reliefrooftop.f90
		rooftop.f90 
		sensorinit.f90 
		sidewall.f90 
		sor3d.f90 
        sort.f90 
		street_intersect.f90 
		streetcanyon.f90 
		surface_coords.f90 
		turbulence_model.f90 
		upwind.f90
		utmll.f90 
		wallbc.f90 
		zd_bisect.f90 )

