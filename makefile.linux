# set path of OPEN_ARRAY
export OPEN_ARRAY=/GPFS/cess/huangxing/openarray_cxx_cess/build

# set path of source code
SRCDIR           = ./src/

OBJDIR           = ./lib/
BINDIR           = ./bin/

EXT_LIB = -L${EXT_PATH}/lib64/
JIT_LIB = ${EXT_PATH}/lib64/

EXE              = GOMO
FC	         = mpiifort -O3 -g -DBOOST_LOG_DYN_LINK -w
FLINKER          = mpiifort -O3 -g -w
CFLAGS	         =
FFLAGS	         =-Wno-tabs -I ${EXT_PATH}/include \
		  -J ${OBJDIR} -I ${OPEN_ARRAY} -g \
		  -fbacktrace -ffree-line-length-0 
CPPFLAGS         =
FPPFLAGS         =
CLEANFILES       = GOMO *.o *.mod *.nc
NP               = 1

OBJ	= ${addprefix ${OBJDIR}/, config.o \
		variables.o dens.o read_init.o  \
		init_fields.o update_initial.o \
		bottom_friction.o get_time.o \
		read_var.o surface_forcing.o \
		lateral_bc.o advct.o baropg.o \
		lateral_viscosity.o advave.o \
		mode_interaction.o bcond1.o \
		external_el.o bcond2_ua.o \
		external_ua.o external_va.o \
		bcond2_va.o external_update.o \
		adjust_uv.o internal_w.o\
		internal_q.o internal_update.o \
		bcond6.o smoth_update.o \
		advt2.o bcond4.o internal_t.o \
		smol_adif.o bcond3_u.o bcond3_v.o \
		internal_u.o internal_v.o print_section.o \
		adjust_ufvf.o gomo.o}

OBJMAIN	= ${OBJ}

${OBJDIR}/%.o: ${SRCDIR}/%.F90
	-${FC} ${FFLAGS} -c -o $@ $<

%.o: %.mod

${EXE}: ${OBJ}
	-${FLINKER} -o ${BINDIR}/${EXE} ${OBJ} \
	${EXT_LIB} -I ${OPEN_ARRAY} -L${OPEN_ARRAY} \
	-lopenarray -lm -ldl -lstdc++ \
	-lboost_program_options -lboost_system -lboost_log -lboost_log_setup -lboost_thread -ljit -lpnetcdf 
clean ::
	-rm lib/*.o
	-rm bin/${EXE}

check :
	-echo ${OBJ}
	-echo ${DM}

