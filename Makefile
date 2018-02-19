# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.7

# Default target executed when no arguments are given to make.
default_target: all

.PHONY : default_target

# Allow only one "make -f Makefile2" at a time, but pass parallelism.
.NOTPARALLEL:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/marlon/Documents/Uni/MasterCHemie/projektmodul/MCLattice2

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/marlon/Documents/Uni/MasterCHemie/projektmodul/MCLattice2

#=============================================================================
# Targets provided globally by CMake.

# Special rule for the target install/strip
install/strip: preinstall
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Installing the project stripped..."
	/usr/bin/cmake -DCMAKE_INSTALL_DO_STRIP=1 -P cmake_install.cmake
.PHONY : install/strip

# Special rule for the target install/strip
install/strip/fast: install/strip

.PHONY : install/strip/fast

# Special rule for the target install/local
install/local: preinstall
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Installing only the local directory..."
	/usr/bin/cmake -DCMAKE_INSTALL_LOCAL_ONLY=1 -P cmake_install.cmake
.PHONY : install/local

# Special rule for the target install/local
install/local/fast: install/local

.PHONY : install/local/fast

# Special rule for the target install
install: preinstall
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Install the project..."
	/usr/bin/cmake -P cmake_install.cmake
.PHONY : install

# Special rule for the target install
install/fast: preinstall/fast
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Install the project..."
	/usr/bin/cmake -P cmake_install.cmake
.PHONY : install/fast

# Special rule for the target list_install_components
list_install_components:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Available install components are: \"Unspecified\""
.PHONY : list_install_components

# Special rule for the target list_install_components
list_install_components/fast: list_install_components

.PHONY : list_install_components/fast

# Special rule for the target rebuild_cache
rebuild_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Running CMake to regenerate build system..."
	/usr/bin/cmake -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR)
.PHONY : rebuild_cache

# Special rule for the target rebuild_cache
rebuild_cache/fast: rebuild_cache

.PHONY : rebuild_cache/fast

# Special rule for the target edit_cache
edit_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "No interactive CMake dialog available..."
	/usr/bin/cmake -E echo No\ interactive\ CMake\ dialog\ available.
.PHONY : edit_cache

# Special rule for the target edit_cache
edit_cache/fast: edit_cache

.PHONY : edit_cache/fast

# The main all target
all: cmake_check_build_system
	$(CMAKE_COMMAND) -E cmake_progress_start /home/marlon/Documents/Uni/MasterCHemie/projektmodul/MCLattice2/CMakeFiles /home/marlon/Documents/Uni/MasterCHemie/projektmodul/MCLattice2/CMakeFiles/progress.marks
	$(MAKE) -f CMakeFiles/Makefile2 all
	$(CMAKE_COMMAND) -E cmake_progress_start /home/marlon/Documents/Uni/MasterCHemie/projektmodul/MCLattice2/CMakeFiles 0
.PHONY : all

# The main clean target
clean:
	$(MAKE) -f CMakeFiles/Makefile2 clean
.PHONY : clean

# The main clean target
clean/fast: clean

.PHONY : clean/fast

# Prepare targets for installation.
preinstall: all
	$(MAKE) -f CMakeFiles/Makefile2 preinstall
.PHONY : preinstall

# Prepare targets for installation.
preinstall/fast:
	$(MAKE) -f CMakeFiles/Makefile2 preinstall
.PHONY : preinstall/fast

# clear depends
depend:
	$(CMAKE_COMMAND) -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 1
.PHONY : depend

#=============================================================================
# Target rules for targets named mclattice2

# Build rule for target.
mclattice2: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 mclattice2
.PHONY : mclattice2

# fast build rule for target.
mclattice2/fast:
	$(MAKE) -f CMakeFiles/mclattice2.dir/build.make CMakeFiles/mclattice2.dir/build
.PHONY : mclattice2/fast

#=============================================================================
# Target rules for targets named enhance

# Build rule for target.
enhance: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 enhance
.PHONY : enhance

# fast build rule for target.
enhance/fast:
	$(MAKE) -f CMakeFiles/enhance.dir/build.make CMakeFiles/enhance.dir/build
.PHONY : enhance/fast

lib/enhance.o: lib/enhance.cpp.o

.PHONY : lib/enhance.o

# target to build an object file
lib/enhance.cpp.o:
	$(MAKE) -f CMakeFiles/enhance.dir/build.make CMakeFiles/enhance.dir/lib/enhance.cpp.o
.PHONY : lib/enhance.cpp.o

lib/enhance.i: lib/enhance.cpp.i

.PHONY : lib/enhance.i

# target to preprocess a source file
lib/enhance.cpp.i:
	$(MAKE) -f CMakeFiles/enhance.dir/build.make CMakeFiles/enhance.dir/lib/enhance.cpp.i
.PHONY : lib/enhance.cpp.i

lib/enhance.s: lib/enhance.cpp.s

.PHONY : lib/enhance.s

# target to generate assembly for a file
lib/enhance.cpp.s:
	$(MAKE) -f CMakeFiles/enhance.dir/build.make CMakeFiles/enhance.dir/lib/enhance.cpp.s
.PHONY : lib/enhance.cpp.s

src/main.o: src/main.cpp.o

.PHONY : src/main.o

# target to build an object file
src/main.cpp.o:
	$(MAKE) -f CMakeFiles/mclattice2.dir/build.make CMakeFiles/mclattice2.dir/src/main.cpp.o
.PHONY : src/main.cpp.o

src/main.i: src/main.cpp.i

.PHONY : src/main.i

# target to preprocess a source file
src/main.cpp.i:
	$(MAKE) -f CMakeFiles/mclattice2.dir/build.make CMakeFiles/mclattice2.dir/src/main.cpp.i
.PHONY : src/main.cpp.i

src/main.s: src/main.cpp.s

.PHONY : src/main.s

# target to generate assembly for a file
src/main.cpp.s:
	$(MAKE) -f CMakeFiles/mclattice2.dir/build.make CMakeFiles/mclattice2.dir/src/main.cpp.s
.PHONY : src/main.cpp.s

src/system/datafile.o: src/system/datafile.cpp.o

.PHONY : src/system/datafile.o

# target to build an object file
src/system/datafile.cpp.o:
	$(MAKE) -f CMakeFiles/mclattice2.dir/build.make CMakeFiles/mclattice2.dir/src/system/datafile.cpp.o
.PHONY : src/system/datafile.cpp.o

src/system/datafile.i: src/system/datafile.cpp.i

.PHONY : src/system/datafile.i

# target to preprocess a source file
src/system/datafile.cpp.i:
	$(MAKE) -f CMakeFiles/mclattice2.dir/build.make CMakeFiles/mclattice2.dir/src/system/datafile.cpp.i
.PHONY : src/system/datafile.cpp.i

src/system/datafile.s: src/system/datafile.cpp.s

.PHONY : src/system/datafile.s

# target to generate assembly for a file
src/system/datafile.cpp.s:
	$(MAKE) -f CMakeFiles/mclattice2.dir/build.make CMakeFiles/mclattice2.dir/src/system/datafile.cpp.s
.PHONY : src/system/datafile.cpp.s

src/system/inputfile.o: src/system/inputfile.cpp.o

.PHONY : src/system/inputfile.o

# target to build an object file
src/system/inputfile.cpp.o:
	$(MAKE) -f CMakeFiles/mclattice2.dir/build.make CMakeFiles/mclattice2.dir/src/system/inputfile.cpp.o
.PHONY : src/system/inputfile.cpp.o

src/system/inputfile.i: src/system/inputfile.cpp.i

.PHONY : src/system/inputfile.i

# target to preprocess a source file
src/system/inputfile.cpp.i:
	$(MAKE) -f CMakeFiles/mclattice2.dir/build.make CMakeFiles/mclattice2.dir/src/system/inputfile.cpp.i
.PHONY : src/system/inputfile.cpp.i

src/system/inputfile.s: src/system/inputfile.cpp.s

.PHONY : src/system/inputfile.s

# target to generate assembly for a file
src/system/inputfile.cpp.s:
	$(MAKE) -f CMakeFiles/mclattice2.dir/build.make CMakeFiles/mclattice2.dir/src/system/inputfile.cpp.s
.PHONY : src/system/inputfile.cpp.s

src/system/lipid.o: src/system/lipid.cpp.o

.PHONY : src/system/lipid.o

# target to build an object file
src/system/lipid.cpp.o:
	$(MAKE) -f CMakeFiles/mclattice2.dir/build.make CMakeFiles/mclattice2.dir/src/system/lipid.cpp.o
.PHONY : src/system/lipid.cpp.o

src/system/lipid.i: src/system/lipid.cpp.i

.PHONY : src/system/lipid.i

# target to preprocess a source file
src/system/lipid.cpp.i:
	$(MAKE) -f CMakeFiles/mclattice2.dir/build.make CMakeFiles/mclattice2.dir/src/system/lipid.cpp.i
.PHONY : src/system/lipid.cpp.i

src/system/lipid.s: src/system/lipid.cpp.s

.PHONY : src/system/lipid.s

# target to generate assembly for a file
src/system/lipid.cpp.s:
	$(MAKE) -f CMakeFiles/mclattice2.dir/build.make CMakeFiles/mclattice2.dir/src/system/lipid.cpp.s
.PHONY : src/system/lipid.cpp.s

src/system/lipidproperties.o: src/system/lipidproperties.cpp.o

.PHONY : src/system/lipidproperties.o

# target to build an object file
src/system/lipidproperties.cpp.o:
	$(MAKE) -f CMakeFiles/mclattice2.dir/build.make CMakeFiles/mclattice2.dir/src/system/lipidproperties.cpp.o
.PHONY : src/system/lipidproperties.cpp.o

src/system/lipidproperties.i: src/system/lipidproperties.cpp.i

.PHONY : src/system/lipidproperties.i

# target to preprocess a source file
src/system/lipidproperties.cpp.i:
	$(MAKE) -f CMakeFiles/mclattice2.dir/build.make CMakeFiles/mclattice2.dir/src/system/lipidproperties.cpp.i
.PHONY : src/system/lipidproperties.cpp.i

src/system/lipidproperties.s: src/system/lipidproperties.cpp.s

.PHONY : src/system/lipidproperties.s

# target to generate assembly for a file
src/system/lipidproperties.cpp.s:
	$(MAKE) -f CMakeFiles/mclattice2.dir/build.make CMakeFiles/mclattice2.dir/src/system/lipidproperties.cpp.s
.PHONY : src/system/lipidproperties.cpp.s

src/system/lipidsystem.o: src/system/lipidsystem.cpp.o

.PHONY : src/system/lipidsystem.o

# target to build an object file
src/system/lipidsystem.cpp.o:
	$(MAKE) -f CMakeFiles/mclattice2.dir/build.make CMakeFiles/mclattice2.dir/src/system/lipidsystem.cpp.o
.PHONY : src/system/lipidsystem.cpp.o

src/system/lipidsystem.i: src/system/lipidsystem.cpp.i

.PHONY : src/system/lipidsystem.i

# target to preprocess a source file
src/system/lipidsystem.cpp.i:
	$(MAKE) -f CMakeFiles/mclattice2.dir/build.make CMakeFiles/mclattice2.dir/src/system/lipidsystem.cpp.i
.PHONY : src/system/lipidsystem.cpp.i

src/system/lipidsystem.s: src/system/lipidsystem.cpp.s

.PHONY : src/system/lipidsystem.s

# target to generate assembly for a file
src/system/lipidsystem.cpp.s:
	$(MAKE) -f CMakeFiles/mclattice2.dir/build.make CMakeFiles/mclattice2.dir/src/system/lipidsystem.cpp.s
.PHONY : src/system/lipidsystem.cpp.s

src/system/mchost.o: src/system/mchost.cpp.o

.PHONY : src/system/mchost.o

# target to build an object file
src/system/mchost.cpp.o:
	$(MAKE) -f CMakeFiles/mclattice2.dir/build.make CMakeFiles/mclattice2.dir/src/system/mchost.cpp.o
.PHONY : src/system/mchost.cpp.o

src/system/mchost.i: src/system/mchost.cpp.i

.PHONY : src/system/mchost.i

# target to preprocess a source file
src/system/mchost.cpp.i:
	$(MAKE) -f CMakeFiles/mclattice2.dir/build.make CMakeFiles/mclattice2.dir/src/system/mchost.cpp.i
.PHONY : src/system/mchost.cpp.i

src/system/mchost.s: src/system/mchost.cpp.s

.PHONY : src/system/mchost.s

# target to generate assembly for a file
src/system/mchost.cpp.s:
	$(MAKE) -f CMakeFiles/mclattice2.dir/build.make CMakeFiles/mclattice2.dir/src/system/mchost.cpp.s
.PHONY : src/system/mchost.cpp.s

src/system/typeproperties.o: src/system/typeproperties.cpp.o

.PHONY : src/system/typeproperties.o

# target to build an object file
src/system/typeproperties.cpp.o:
	$(MAKE) -f CMakeFiles/mclattice2.dir/build.make CMakeFiles/mclattice2.dir/src/system/typeproperties.cpp.o
.PHONY : src/system/typeproperties.cpp.o

src/system/typeproperties.i: src/system/typeproperties.cpp.i

.PHONY : src/system/typeproperties.i

# target to preprocess a source file
src/system/typeproperties.cpp.i:
	$(MAKE) -f CMakeFiles/mclattice2.dir/build.make CMakeFiles/mclattice2.dir/src/system/typeproperties.cpp.i
.PHONY : src/system/typeproperties.cpp.i

src/system/typeproperties.s: src/system/typeproperties.cpp.s

.PHONY : src/system/typeproperties.s

# target to generate assembly for a file
src/system/typeproperties.cpp.s:
	$(MAKE) -f CMakeFiles/mclattice2.dir/build.make CMakeFiles/mclattice2.dir/src/system/typeproperties.cpp.s
.PHONY : src/system/typeproperties.cpp.s

# Help Target
help:
	@echo "The following are some of the valid targets for this Makefile:"
	@echo "... all (the default if no target is provided)"
	@echo "... clean"
	@echo "... depend"
	@echo "... install/strip"
	@echo "... install/local"
	@echo "... install"
	@echo "... list_install_components"
	@echo "... mclattice2"
	@echo "... rebuild_cache"
	@echo "... enhance"
	@echo "... edit_cache"
	@echo "... lib/enhance.o"
	@echo "... lib/enhance.i"
	@echo "... lib/enhance.s"
	@echo "... src/main.o"
	@echo "... src/main.i"
	@echo "... src/main.s"
	@echo "... src/system/datafile.o"
	@echo "... src/system/datafile.i"
	@echo "... src/system/datafile.s"
	@echo "... src/system/inputfile.o"
	@echo "... src/system/inputfile.i"
	@echo "... src/system/inputfile.s"
	@echo "... src/system/lipid.o"
	@echo "... src/system/lipid.i"
	@echo "... src/system/lipid.s"
	@echo "... src/system/lipidproperties.o"
	@echo "... src/system/lipidproperties.i"
	@echo "... src/system/lipidproperties.s"
	@echo "... src/system/lipidsystem.o"
	@echo "... src/system/lipidsystem.i"
	@echo "... src/system/lipidsystem.s"
	@echo "... src/system/mchost.o"
	@echo "... src/system/mchost.i"
	@echo "... src/system/mchost.s"
	@echo "... src/system/typeproperties.o"
	@echo "... src/system/typeproperties.i"
	@echo "... src/system/typeproperties.s"
.PHONY : help



#=============================================================================
# Special targets to cleanup operation of make.

# Special rule to run CMake to check the build system integrity.
# No rule that depends on this can have commands that come from listfiles
# because they might be regenerated.
cmake_check_build_system:
	$(CMAKE_COMMAND) -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 0
.PHONY : cmake_check_build_system

