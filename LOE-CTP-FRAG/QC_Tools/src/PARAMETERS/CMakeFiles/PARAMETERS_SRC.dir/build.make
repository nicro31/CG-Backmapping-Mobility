# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 2.8

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

# The program to use to edit the cache.
CMAKE_EDIT_COMMAND = /usr/bin/ccmake

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/x_nicro/LOE-CTP-FRAG/QC_Tools

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/x_nicro/LOE-CTP-FRAG/QC_Tools

# Include any dependencies generated for this target.
include src/PARAMETERS/CMakeFiles/PARAMETERS_SRC.dir/depend.make

# Include the progress variables for this target.
include src/PARAMETERS/CMakeFiles/PARAMETERS_SRC.dir/progress.make

# Include the compile flags for this target's objects.
include src/PARAMETERS/CMakeFiles/PARAMETERS_SRC.dir/flags.make

src/PARAMETERS/CMakeFiles/PARAMETERS_SRC.dir/parameters.cpp.o: src/PARAMETERS/CMakeFiles/PARAMETERS_SRC.dir/flags.make
src/PARAMETERS/CMakeFiles/PARAMETERS_SRC.dir/parameters.cpp.o: src/PARAMETERS/parameters.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/x_nicro/LOE-CTP-FRAG/QC_Tools/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object src/PARAMETERS/CMakeFiles/PARAMETERS_SRC.dir/parameters.cpp.o"
	cd /home/x_nicro/LOE-CTP-FRAG/QC_Tools/src/PARAMETERS && /software/sse/wrappers/icpc   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/PARAMETERS_SRC.dir/parameters.cpp.o -c /home/x_nicro/LOE-CTP-FRAG/QC_Tools/src/PARAMETERS/parameters.cpp

src/PARAMETERS/CMakeFiles/PARAMETERS_SRC.dir/parameters.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/PARAMETERS_SRC.dir/parameters.cpp.i"
	cd /home/x_nicro/LOE-CTP-FRAG/QC_Tools/src/PARAMETERS && /software/sse/wrappers/icpc  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/x_nicro/LOE-CTP-FRAG/QC_Tools/src/PARAMETERS/parameters.cpp > CMakeFiles/PARAMETERS_SRC.dir/parameters.cpp.i

src/PARAMETERS/CMakeFiles/PARAMETERS_SRC.dir/parameters.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/PARAMETERS_SRC.dir/parameters.cpp.s"
	cd /home/x_nicro/LOE-CTP-FRAG/QC_Tools/src/PARAMETERS && /software/sse/wrappers/icpc  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/x_nicro/LOE-CTP-FRAG/QC_Tools/src/PARAMETERS/parameters.cpp -o CMakeFiles/PARAMETERS_SRC.dir/parameters.cpp.s

src/PARAMETERS/CMakeFiles/PARAMETERS_SRC.dir/parameters.cpp.o.requires:
.PHONY : src/PARAMETERS/CMakeFiles/PARAMETERS_SRC.dir/parameters.cpp.o.requires

src/PARAMETERS/CMakeFiles/PARAMETERS_SRC.dir/parameters.cpp.o.provides: src/PARAMETERS/CMakeFiles/PARAMETERS_SRC.dir/parameters.cpp.o.requires
	$(MAKE) -f src/PARAMETERS/CMakeFiles/PARAMETERS_SRC.dir/build.make src/PARAMETERS/CMakeFiles/PARAMETERS_SRC.dir/parameters.cpp.o.provides.build
.PHONY : src/PARAMETERS/CMakeFiles/PARAMETERS_SRC.dir/parameters.cpp.o.provides

src/PARAMETERS/CMakeFiles/PARAMETERS_SRC.dir/parameters.cpp.o.provides.build: src/PARAMETERS/CMakeFiles/PARAMETERS_SRC.dir/parameters.cpp.o

# Object files for target PARAMETERS_SRC
PARAMETERS_SRC_OBJECTS = \
"CMakeFiles/PARAMETERS_SRC.dir/parameters.cpp.o"

# External object files for target PARAMETERS_SRC
PARAMETERS_SRC_EXTERNAL_OBJECTS =

src/PARAMETERS/libPARAMETERS_SRC.a: src/PARAMETERS/CMakeFiles/PARAMETERS_SRC.dir/parameters.cpp.o
src/PARAMETERS/libPARAMETERS_SRC.a: src/PARAMETERS/CMakeFiles/PARAMETERS_SRC.dir/build.make
src/PARAMETERS/libPARAMETERS_SRC.a: src/PARAMETERS/CMakeFiles/PARAMETERS_SRC.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX static library libPARAMETERS_SRC.a"
	cd /home/x_nicro/LOE-CTP-FRAG/QC_Tools/src/PARAMETERS && $(CMAKE_COMMAND) -P CMakeFiles/PARAMETERS_SRC.dir/cmake_clean_target.cmake
	cd /home/x_nicro/LOE-CTP-FRAG/QC_Tools/src/PARAMETERS && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/PARAMETERS_SRC.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/PARAMETERS/CMakeFiles/PARAMETERS_SRC.dir/build: src/PARAMETERS/libPARAMETERS_SRC.a
.PHONY : src/PARAMETERS/CMakeFiles/PARAMETERS_SRC.dir/build

src/PARAMETERS/CMakeFiles/PARAMETERS_SRC.dir/requires: src/PARAMETERS/CMakeFiles/PARAMETERS_SRC.dir/parameters.cpp.o.requires
.PHONY : src/PARAMETERS/CMakeFiles/PARAMETERS_SRC.dir/requires

src/PARAMETERS/CMakeFiles/PARAMETERS_SRC.dir/clean:
	cd /home/x_nicro/LOE-CTP-FRAG/QC_Tools/src/PARAMETERS && $(CMAKE_COMMAND) -P CMakeFiles/PARAMETERS_SRC.dir/cmake_clean.cmake
.PHONY : src/PARAMETERS/CMakeFiles/PARAMETERS_SRC.dir/clean

src/PARAMETERS/CMakeFiles/PARAMETERS_SRC.dir/depend:
	cd /home/x_nicro/LOE-CTP-FRAG/QC_Tools && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/x_nicro/LOE-CTP-FRAG/QC_Tools /home/x_nicro/LOE-CTP-FRAG/QC_Tools/src/PARAMETERS /home/x_nicro/LOE-CTP-FRAG/QC_Tools /home/x_nicro/LOE-CTP-FRAG/QC_Tools/src/PARAMETERS /home/x_nicro/LOE-CTP-FRAG/QC_Tools/src/PARAMETERS/CMakeFiles/PARAMETERS_SRC.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/PARAMETERS/CMakeFiles/PARAMETERS_SRC.dir/depend
