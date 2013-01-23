#!/bin/sh

built_sources=""

headers=`find .. -name "*.h" -type f | sort`

for header_with_path in $headers ; do
    
    #echo $header_with_path
    header=`basename $header_with_path`
    #echo $header
    built_sources="$built_sources $header"
done

specializations=`find .. -name "*specializations" -type f | sort`

for specialization_with_path in $specializations ; do
    
    #echo $specialization_with_path
    specialization=`basename $specialization_with_path`
    #echo $specialization
    built_sources="$built_sources $specialization"
done

cat <<EOF > Makefile.am
# Note - this file is automatically generated by $0 
# do not edit manually

#
# include the magic script!
EXTRA_DIST = rebuild_makefile.sh

BUILT_SOURCES = $built_sources

DISTCLEANFILES = \$(BUILT_SOURCES)

EOF



# handle contrib directly
cat <<EOF >> Makefile.am
#
# contrib rules
if LIBMESH_ENABLE_FPARSER

fparser.hh: \$(top_srcdir)/contrib/fparser/fparser.hh
	\$(AM_V_GEN)\$(LN_S) \$(top_srcdir)/contrib/fparser/fparser.hh fparser.hh

  BUILT_SOURCES  += fparser.hh
  DISTCLEANFILES += fparser.hh

endif

if LIBMESH_ENABLE_NANOFLANN

nanoflann.hpp: \$(top_srcdir)/contrib/nanoflann/include/nanoflann.hpp
	\$(AM_V_GEN)\$(LN_S) \$(top_srcdir)/contrib/nanoflann/include/nanoflann.hpp nanoflann.hpp

  BUILT_SOURCES  += nanoflann.hpp
  DISTCLEANFILES += nanoflann.hpp

endif

if LIBMESH_ENABLE_EXODUS

exodusII.h: \$(top_srcdir)/contrib/exodusii/Lib/include/exodusII.h
	\$(AM_V_GEN)\$(LN_S) \$(top_srcdir)/contrib/exodusii/Lib/include/exodusII.h exodusII.h

  BUILT_SOURCES  += exodusII.h
  DISTCLEANFILES += exodusII.h

exodusII_ext.h: \$(top_srcdir)/contrib/exodusii/Lib/include/exodusII_ext.h
	\$(AM_V_GEN)\$(LN_S) \$(top_srcdir)/contrib/exodusii/Lib/include/exodusII_ext.h exodusII_ext.h

  BUILT_SOURCES  += exodusII_ext.h
  DISTCLEANFILES += exodusII_ext.h

netcdf.h: \$(top_srcdir)/contrib/netcdf/Lib/netcdf.h
	\$(AM_V_GEN)\$(LN_S) \$(top_srcdir)/contrib/netcdf/Lib/netcdf.h netcdf.h

  BUILT_SOURCES  += netcdf.h
  DISTCLEANFILES += netcdf.h

endif

EOF



# handle libmesh_config.h
cat <<EOF >> Makefile.am
#
# libmesh_config.h rule
libmesh_config.h: \$(top_builddir)/include/libmesh_config.h
	\$(AM_V_GEN)\$(LN_S) \$(top_builddir)/include/libmesh_config.h libmesh_config.h

  BUILT_SOURCES  += libmesh_config.h
  DISTCLEANFILES += libmesh_config.h

EOF



# now automatically handle our headers
cat <<EOF >> Makefile.am
#
# libMesh header rules
EOF
for header_with_path in $headers $specializations ; do  
    header=`basename $header_with_path`
    source=`echo $header_with_path | sed 's/../$(top_srcdir)\/include/' -`
    #echo $source
    cat <<EOF >> Makefile.am
$header: $source
	\$(AM_V_GEN)\$(LN_S) $source $header

EOF
done
#cat Makefile.am
