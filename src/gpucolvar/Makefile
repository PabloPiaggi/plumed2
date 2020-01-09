USE=core reference tools

#generic makefile
#
# This is a generic makefile which will work for any
# subpackage which just depends on config/
# It just compiles objects locally
# 

# include the machine dependent configuration
ifneq ($(MAKECMDGOALS),clean)
  -include ../../Makefile.conf
endif

# if machine dependent configuration has been found:
ifdef GCCDEP

# source for plumed kernel
SRC=$(wildcard *.cpp)
SRCCU=$(wildcard *.cu)

GPUCOLVAR_LIB_NAME=libgpucolvar.so

# corresponding objects:
OBJ_CPP=$(SRC:.cpp=.o)
OBJ_CUDA=$(SRCCU:.cu=.o)
OBJ_CUDA_DLINK=$(SRCCU:.cu=.dlink.o)
OBJ=$(OBJ_CPP)
OBJ+=$(OBJ_CUDA)
OBJ+=$(OBJ_CUDA_DLINK)
OBJ+=$(GPUCOLVAR_LIB_NAME)

# dependency files for all c++ sources
DEP=$(addprefix deps/,$(SRC:.cpp=.d))
DEP+=$(addprefix deps/,$(SRC:.cu=.d))

.PHONY: all dirslinks obj lib install

# default target:
# update all links then build this module
all: dirslinks
	@echo Preparing for make
	$(PREMAKE)
	@echo "*** Compile local objects ***"
	$(MAKE) obj


# build links for all modules
dirslinks:
	make -C ../lib dirslinks

# buils objects in this module
obj: $(OBJ)

# build everything
lib:
	cd ../lib ; make

.PHONY: install
install:
	cd ../lib ; make install

ifneq ($(disable_dependency_tracking),yes)
-include $(DEP)
endif

# if machine dependent configuration has not been found:
else

.PHONY: error
error:
	@echo No configuration available
	@echo First run ./configure in the root directory
endif

# this target is available anyway

.PHONY: clean links
clean:
	@echo Clean $(CURDIR)
	../maketools/cleanlinks.sh
	rm -fr deps links
	rm -f $(CLEANLIST)

# build links for this module
links: Makefile
	@echo Verify modules and create links in $(CURDIR)
	@for dir in $(USE) ; do ../maketools/check_module.sh $$dir ; ../maketools/makelinks.sh $$dir ;  done

.PHONY: cppcheck
cppcheck:
	../maketools/codecheck --cppcheck --local

.PHONY: plumedcheck
plumedcheck:
	../maketools/codecheck --plumedcheck --local

.PHONY: codecheck
codecheck:
	../maketools/codecheck --local

.PHONY: astyle
astyle:
	cd ../ ; ./astyle.sh $(CURDIR)

.PHONY: show_used_modules
show_used_modules:
	@echo used_modules $(USE)



# generic makefile rules
# this compiler is used for dependencies generation
# if unset, use normal CXX
ifndef CXXDEP
CXXDEP=$(CXX)
CCDEP=$(CC)
endif

ifndef VERBOSE
AT=@
endif

# rule for c++
%.o: %.cpp
	@echo "compiling " $*.cpp
	@-test -d deps || mkdir deps
ifneq ($(disable_dependency_tracking),yes)
	$(AT)$(CXXDEP) -c -MM -MF$*.d $(CPPFLAGS) $(ADDCPPFLAGS) $(CXXFLAGS) $*.cpp
	@cp -f $*.d $*.d.tmp
	@sed -e 's/.*://' -e 's/\\$$//' < $*.d.tmp | fmt -1 | \
	 sed -e 's/^ *//' -e 's/$$/:/' >> $*.d
	@mv $*.d deps/$*.d
	@rm -f $*.d.tmp
endif
ifndef XLF
	$(AT)$(CXX) -c $(CPPFLAGS) $(ADDCPPFLAGS) $(CXXFLAGS) $*.cpp -o $*.o
endif

# rule for c
%.o: %.c
	@echo "compiling " $*.c
	@-test -d deps || mkdir deps
ifneq ($(disable_dependency_tracking),yes)
	$(AT)$(CXXDEP) -c -MM -MF$*.d $(CPPFLAGS) $(ADDCPPFLAGS) $(CFLAGS) $*.c -o $*.o
	@cp -f $*.d $*.d.tmp
	@sed -e 's/.*://' -e 's/\\$$//' < $*.d.tmp | fmt -1 | \
	 sed -e 's/^ *//' -e 's/$$/:/' >> $*.d
	@mv $*.d deps/$*.d
	@rm -f $*.d.tmp
endif
ifndef XLF
	$(AT)$(CC) -c $(CPPFLAGS) $(ADDCPPFLAGS) $(CFLAGS) $*.c -o $*.o
endif

ifndef CUDA
CUDA=nvcc
endif

ifndef CUDAFLAGS
CUDAFLAGS=-arch=sm_60 -std=c++11 -shared -Xcompiler -fPIC -g -G
endif

# rule for cuda
%.o: %.cu
	@echo "compiling " $*.cu
	@-test -d deps || mkdir deps
ifneq ($(disable_dependency_tracking),yes)
	$(AT)$(CUDA) -c -MD $(CUDAFLAGS) -dc $*.cu -o $*.o -MF $*.d 
	@cp -f $*.d $*.d.tmp
	@sed -e 's/.*://' -e 's/\\$$//' < $*.d.tmp | fmt -1 | \
	 sed -e 's/^ *//' -e 's/$$/:/' >> $*.d
	@mv $*.d deps/$*.d
	@rm -f $*.d.tmp
endif
ifndef XLF
	$(AT)$(CUDA) $(CUDAFLAGS) -c $*.cu -o $*.o
endif

# Device link step
%.dlink.o: %.o
	@echo "Device link step " $*.o
	$(AT)$(CUDA) $(CUDAFLAGS) -dlink -o $*.dlink.o $*.o

# Create library
$(GPUCOLVAR_LIB_NAME): $(OBJ_CUDA) $(OBJ_CUDA_DLINK)
	@echo "Creating library " $@
	$(AT)$(CUDA) $(CUDAFLAGS) -o $@ $^

#
CLEANLIST=*~ *.tmp *.xxd *.o *.so *.dylib *.gcda *.gcno *.a

