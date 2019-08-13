AUTOMAKE_OPTIONS = foreign

oplevel0 =  -O3 -g -w
FCFLAGS += -fbacktrace -ffree-line-length-0 -I/usr/local/include -I /usr/include $(oplevel0)
CFLAGS += -I/usr/local/include -I /usr/include $(oplevel0)
CXXFLAGS += -I/usr/local/include -I /usr/include $(oplevel0)

bin_PROGRAMS = $(top_builddir)/bin/GOMO
__top_builddir__bin_GOMO_SOURCES = src/config.F90
__top_builddir__bin_GOMO_SOURCES += src/variables.F90
__top_builddir__bin_GOMO_SOURCES += src/dens.F90
__top_builddir__bin_GOMO_SOURCES += src/read_init.F90
__top_builddir__bin_GOMO_SOURCES += src/init_fields.F90
__top_builddir__bin_GOMO_SOURCES += src/update_initial.F90
__top_builddir__bin_GOMO_SOURCES += src/bottom_friction.F90
__top_builddir__bin_GOMO_SOURCES += src/get_time.F90
__top_builddir__bin_GOMO_SOURCES += src/read_var.F90
__top_builddir__bin_GOMO_SOURCES += src/surface_forcing.F90
__top_builddir__bin_GOMO_SOURCES += src/lateral_bc.F90
__top_builddir__bin_GOMO_SOURCES += src/advct.F90
__top_builddir__bin_GOMO_SOURCES += src/baropg.F90
__top_builddir__bin_GOMO_SOURCES += src/lateral_viscosity.F90
__top_builddir__bin_GOMO_SOURCES += src/advave.F90
__top_builddir__bin_GOMO_SOURCES += src/mode_interaction.F90
__top_builddir__bin_GOMO_SOURCES += src/bcond1.F90
__top_builddir__bin_GOMO_SOURCES += src/external_el.F90
__top_builddir__bin_GOMO_SOURCES += src/bcond2_ua.F90
__top_builddir__bin_GOMO_SOURCES += src/external_ua.F90
__top_builddir__bin_GOMO_SOURCES += src/external_va.F90
__top_builddir__bin_GOMO_SOURCES += src/bcond2_va.F90
__top_builddir__bin_GOMO_SOURCES += src/external_update.F90
__top_builddir__bin_GOMO_SOURCES += src/adjust_uv.F90
__top_builddir__bin_GOMO_SOURCES += src/internal_w.F90
__top_builddir__bin_GOMO_SOURCES += src/internal_q.F90
__top_builddir__bin_GOMO_SOURCES += src/internal_update.F90
__top_builddir__bin_GOMO_SOURCES += src/bcond6.F90
__top_builddir__bin_GOMO_SOURCES += src/smoth_update.F90
__top_builddir__bin_GOMO_SOURCES += src/advt2.F90
__top_builddir__bin_GOMO_SOURCES += src/bcond4.F90
__top_builddir__bin_GOMO_SOURCES += src/internal_t.F90
__top_builddir__bin_GOMO_SOURCES += src/smol_adif.F90
__top_builddir__bin_GOMO_SOURCES += src/bcond3_u.F90
__top_builddir__bin_GOMO_SOURCES += src/bcond3_v.F90
__top_builddir__bin_GOMO_SOURCES += src/internal_u.F90
__top_builddir__bin_GOMO_SOURCES += src/internal_v.F90
__top_builddir__bin_GOMO_SOURCES += src/print_section.F90
__top_builddir__bin_GOMO_SOURCES += src/adjust_ufvf.F90
__top_builddir__bin_GOMO_SOURCES += src/gomo.F90

__top_builddir__bin_GOMO_LDADD = -lopenarray -lm -ldl -lstdc++ -L. -lstdc++ -lgfortran -lm -ldl -lpthread ${AM_LDFLAGS}

src/adjust_ufvf.o : src/adjust_ufvf.F90 src/variables.o src/config.o
src/adjust_uv.o : src/adjust_uv.F90 src/config.o src/variables.o
src/advave.o : src/advave.F90 src/config.o src/variables.o
src/advct.o : src/advct.F90 src/variables.o src/config.o
src/advt2.o : src/advt2.F90 src/config.o src/variables.o
src/baropg.o : src/baropg.F90 src/config.o src/variables.o
src/bcond1.o : src/bcond1.F90 src/config.o src/variables.o
src/bcond2_ua.o : src/bcond2_ua.F90 src/config.o src/variables.o
src/bcond2_va.o : src/bcond2_va.F90 src/config.o src/variables.o
src/bcond3_u.o : src/bcond3_u.F90 src/config.o src/variables.o
src/bcond3_v.o : src/bcond3_v.F90 src/config.o src/variables.o
src/bcond4.o : src/bcond4.F90 src/variables.o src/config.o
src/bcond6.o : src/bcond6.F90 src/config.o src/variables.o
src/bottom_friction.o : src/bottom_friction.F90 src/config.o src/variables.o
src/config.o : src/config.F90
src/constants.o : src/constants.F90
src/dens.o : src/dens.F90 src/config.o src/variables.o
src/external_el.o : src/external_el.F90 src/config.o src/variables.o
src/external_ua.o : src/external_ua.F90 src/config.o src/variables.o
src/external_update.o : src/external_update.F90 src/variables.o src/config.o
src/external_va.o : src/external_va.F90 src/config.o src/variables.o
src/get_time.o : src/get_time.F90 src/config.o src/variables.o
src/gomo.o : src/gomo.F90 src/variables.o src/config.o
src/init_fields.o : src/init_fields.F90 src/variables.o
src/internal_q.o : src/internal_q.F90 src/variables.o src/config.o
src/internal_t.o : src/internal_t.F90 src/config.o src/variables.o
src/internal_u.o : src/internal_u.F90 src/variables.o src/config.o
src/internal_update.o : src/internal_update.F90 src/variables.o src/config.o
src/internal_v.o : src/internal_v.F90 src/variables.o src/config.o
src/internal_w.o : src/internal_w.F90 src/variables.o src/config.o
src/lateral_bc.o : src/lateral_bc.F90 src/config.o src/variables.o
src/lateral_viscosity.o : src/lateral_viscosity.F90 src/variables.o src/config.o
src/mode_interaction.o : src/mode_interaction.F90 src/variables.o src/config.o
src/print_section.o : src/print_section.F90 src/variables.o src/config.o
src/read_init.o : src/read_init.F90 src/config.o src/variables.o
src/read_var.o : src/read_var.F90 src/config.o src/variables.o
src/read_wind.o : src/read_wind.F90 src/config.o src/variables.o
src/smol_adif.o : src/smol_adif.F90 src/variables.o src/config.o
src/smoth_update.o : src/smoth_update.F90 src/config.o src/variables.o
src/surface_forcing.o : src/surface_forcing.F90 src/config.o src/variables.o
src/update_initial.o : src/update_initial.F90 src/config.o src/variables.o
src/variables.o : src/variables.F90 src/config.o

src:
	git archive --format=tar.gz --prefix=${PACKAGE_NAME}-${PACKAGE_VERSION}/ v${PACKAGE_VERSION} > ${PACKAGE_NAME}-${PACKAGE_VERSION}.tar.gz

bin: GOMO
	zip ${PACKAGE_NAME}-bin-release-${PACKAGE_VERSION}.zip $^
