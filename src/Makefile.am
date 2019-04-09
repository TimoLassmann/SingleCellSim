# Build the library in the hand subdirectory first.


LIB_TLDEVELDIR = @LIB_TLDEVELDIR@

MYLIBDIRS = ../tldevel/libthrpool.a \
../tldevel/libtldevel.a


MYINCDIRS = -I${top_builddir}/${LIB_TLDEVELDIR} \
-I${top_srcdir}/${LIB_TLDEVELDIR}







AM_CPPFLAGS=${MYINCDIRS}

bin_PROGRAMS = scs
scs_SOURCES = \
scs.h \
user_interface.c \
user_interface.h \
singleCellSim.c \
gigp.h \
gigp.c \
matrix_io.c \
outdir.h \
outdir.c

scs_LDADD = ${MYLIBDIRS}

