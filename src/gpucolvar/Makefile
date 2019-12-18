USE=core reference tools

# generic makefile
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

# corresponding objects:
OBJ=$(SRC:.cpp=.o)
OBJ+=$(SRCCU:.cu=.o)

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
ifndef CUDA
CUDA=nvcc
endif

ifndef CUDAFLAGS
CUDAFLAGS=-arch=sm_60 -std=c++11
endif

ifndef VERBOSE
AT=@
endif

# rule for cuda
%.o: %.cu
	@echo "compiling " $*.cu
	@-test -d deps || mkdir deps
	echo "$(AT)$(CUDA) $(CUDAFLAGS) -c $*.cu -o $*.o"
	$(AT)$(CUDA) $(CUDAFLAGS) -c $*.cu -o $*.o


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

#
CLEANLIST=*~ *.tmp *.xxd *.o *.so *.dylib *.gcda *.gcno *.a

