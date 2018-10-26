# You will likely have to change the following two lines
module_LTLIBRARIES = skew_lattice.la
skew_lattice_la_SOURCES = skew_lattice.c

# The rest is quite generic unless your module uses extra libraries
ACLOCAL_AMFLAGS = -I m4
moduledir = @GWYDDION_MODULE_DIR@
AM_CPPFLAGS = -I$(top_srcdir) -DG_LOG_DOMAIN=\"Module\" @GWYDDION_CFLAGS@
AM_CFLAGS = @WARNING_CFLAGS@ @HOST_CFLAGS@
AM_LDFLAGS = -avoid-version -module @HOST_LDFLAGS@ @GWYDDION_LIBS@
