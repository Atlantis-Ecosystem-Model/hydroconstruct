AUTOMAKE_OPTIONS = foreign

OPT = -msse2 -mfpmath=sse

CC = gcc
DBG = -g
DEFINES = -D PROJ4 -D SMALL_LOG_FILE_SIZE

# You will need to modify the pathway to sjwlib as appropaiet - it does not seem to be being picked up from the pass through from "pkg-config --cflags atlantisMerged"
# Use this one from Beth's machine as an example, you basically need to give it the relative path
# from your Hydroconstruct to your sjwlib in your Atlantis code directly
INCLUDES = $(shell pkg-config --cflags atlantisMerged)\
    -I$(top_srcdir)/../../AtlantisCurrentTrunk/atlantis/sjwlib/include

HydroConstruct_LDADD   =  -lsjwlib -lm

#LIBDIRS = -L$(top_srcdir)/../../AtlantisCurrentTrunk/atlantis/sjwlib -I.
#HydroConstruct_LDFLAGS = $(LIBDIRS) $(HydroConstruct_DEPENDENCIES)
#HydroConstruct_DEPENDENCIES =   $(top_srcdir)/../Atlantis/trunk/atlantis/sjwlib/libsjwlib.a  

WARN =
WARN += -Wall -Wextra
WARN += -Wno-pointer-arith
WARN += -Wstrict-prototypes
WARN += -Werror 
WARN += -Wunused-variable
WARN += -Wuninitialized
WARN += -Wformat 


AM_CFLAGS = $(DEFINES) $(INCLUDE) $(DBG) $(OPT) $(EF) $(PG) $(WARN)

bin_PROGRAMS   = HydroConstruct
HydroConstruct_SOURCES = calcexchange.c getaltinputdata.c getcdfinputdata.c getcdftsinputdata.c getflatinputdata.c getinputdata.c \
	getROMScdfinputdata.c getROMScdftsinputdata.c mainexchange.c padtimeseries.c readnetcdf.c writeoutput.c \
	get2in1ROMScdfinputdata.c getEccalFlatFormatData.c checkTransport.c unixUtils.c getGuamData.c


