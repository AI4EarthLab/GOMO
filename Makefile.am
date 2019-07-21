AUTOMAKE_OPTIONS = foreign

oplevel0 =  -O3 -g -w
FCFLAGS += -I/usr/local/include -I /usr/include $(oplevel0)
CFLAGS += -I/usr/local/include -I /usr/include $(oplevel0)
CXXFLAGS += -I/usr/local/include -I /usr/include $(oplevel0)

bin_PROGRAMS = gomo
gomo_SOURCES = src/config.F90
gomo_SOURCES += src/variables.F90
gomo_SOURCES += src/dens.F90
gomo_SOURCES += src/read_init.F90
gomo_SOURCES += src/init_fields.F90
gomo_SOURCES += src/update_initial.F90
gomo_SOURCES += src/bottom_friction.F90
gomo_SOURCES += src/get_time.F90
gomo_SOURCES += src/read_var.F90
gomo_SOURCES += src/surface_forcing.F90
gomo_SOURCES += src/lateral_bc.F90
gomo_SOURCES += src/advct.F90
gomo_SOURCES += src/baropg.F90
gomo_SOURCES += src/lateral_viscosity.F90
gomo_SOURCES += src/advave.F90
gomo_SOURCES += src/mode_interaction.F90
gomo_SOURCES += src/bcond1.F90
gomo_SOURCES += src/external_el.F90
gomo_SOURCES += src/bcond2_ua.F90
gomo_SOURCES += src/external_ua.F90
gomo_SOURCES += src/external_va.F90
gomo_SOURCES += src/bcond2_va.F90
gomo_SOURCES += src/external_update.F90
gomo_SOURCES += src/adjust_uv.F90
gomo_SOURCES += src/internal_w.F90
gomo_SOURCES += src/internal_q.F90
gomo_SOURCES += src/internal_update.F90
gomo_SOURCES += src/bcond6.F90
gomo_SOURCES += src/smoth_update.oF90
gomo_SOURCES += src/advt2.o bcond4.F90
gomo_SOURCES += src/internal_t.F90
gomo_SOURCES += src/smol_adif.F90
gomo_SOURCES += src/bcond3_u.F90
gomo_SOURCES += src/bcond3_v.F90
gomo_SOURCES += src/internal_u.F90
gomo_SOURCES += src/internal_v.F90
gomo_SOURCES += src/print_section.F90
gomo_SOURCES += src/adjust_ufvf.F90
gomo_SOURCES += src/gomo.F90

gomo_LDADD = -lgomo -lm -ldl -lstdc++ -L. -lstdc++ -lgfortran -lm -ldl -lpthread ${AM_LDFLAGS}

src:
	git archive --format=tar.gz --prefix=${PACKAGE_NAME}-${PACKAGE_VERSION}/ v${PACKAGE_VERSION} > ${PACKAGE_NAME}-${PACKAGE_VERSION}.tar.gz

bin: gomo
	zip ${PACKAGE_NAME}-bin-release-${PACKAGE_VERSION}.zip $^
