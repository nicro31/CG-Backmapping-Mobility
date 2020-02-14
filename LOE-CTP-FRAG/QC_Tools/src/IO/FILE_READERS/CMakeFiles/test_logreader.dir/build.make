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
include src/IO/FILE_READERS/CMakeFiles/test_logreader.dir/depend.make

# Include the progress variables for this target.
include src/IO/FILE_READERS/CMakeFiles/test_logreader.dir/progress.make

# Include the compile flags for this target's objects.
include src/IO/FILE_READERS/CMakeFiles/test_logreader.dir/flags.make

src/IO/FILE_READERS/CMakeFiles/test_logreader.dir/test_logreader.cpp.o: src/IO/FILE_READERS/CMakeFiles/test_logreader.dir/flags.make
src/IO/FILE_READERS/CMakeFiles/test_logreader.dir/test_logreader.cpp.o: src/IO/FILE_READERS/test_logreader.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/x_nicro/LOE-CTP-FRAG/QC_Tools/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object src/IO/FILE_READERS/CMakeFiles/test_logreader.dir/test_logreader.cpp.o"
	cd /home/x_nicro/LOE-CTP-FRAG/QC_Tools/src/IO/FILE_READERS && /software/sse/wrappers/icpc   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/test_logreader.dir/test_logreader.cpp.o -c /home/x_nicro/LOE-CTP-FRAG/QC_Tools/src/IO/FILE_READERS/test_logreader.cpp

src/IO/FILE_READERS/CMakeFiles/test_logreader.dir/test_logreader.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/test_logreader.dir/test_logreader.cpp.i"
	cd /home/x_nicro/LOE-CTP-FRAG/QC_Tools/src/IO/FILE_READERS && /software/sse/wrappers/icpc  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/x_nicro/LOE-CTP-FRAG/QC_Tools/src/IO/FILE_READERS/test_logreader.cpp > CMakeFiles/test_logreader.dir/test_logreader.cpp.i

src/IO/FILE_READERS/CMakeFiles/test_logreader.dir/test_logreader.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/test_logreader.dir/test_logreader.cpp.s"
	cd /home/x_nicro/LOE-CTP-FRAG/QC_Tools/src/IO/FILE_READERS && /software/sse/wrappers/icpc  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/x_nicro/LOE-CTP-FRAG/QC_Tools/src/IO/FILE_READERS/test_logreader.cpp -o CMakeFiles/test_logreader.dir/test_logreader.cpp.s

src/IO/FILE_READERS/CMakeFiles/test_logreader.dir/test_logreader.cpp.o.requires:
.PHONY : src/IO/FILE_READERS/CMakeFiles/test_logreader.dir/test_logreader.cpp.o.requires

src/IO/FILE_READERS/CMakeFiles/test_logreader.dir/test_logreader.cpp.o.provides: src/IO/FILE_READERS/CMakeFiles/test_logreader.dir/test_logreader.cpp.o.requires
	$(MAKE) -f src/IO/FILE_READERS/CMakeFiles/test_logreader.dir/build.make src/IO/FILE_READERS/CMakeFiles/test_logreader.dir/test_logreader.cpp.o.provides.build
.PHONY : src/IO/FILE_READERS/CMakeFiles/test_logreader.dir/test_logreader.cpp.o.provides

src/IO/FILE_READERS/CMakeFiles/test_logreader.dir/test_logreader.cpp.o.provides.build: src/IO/FILE_READERS/CMakeFiles/test_logreader.dir/test_logreader.cpp.o

# Object files for target test_logreader
test_logreader_OBJECTS = \
"CMakeFiles/test_logreader.dir/test_logreader.cpp.o"

# External object files for target test_logreader
test_logreader_EXTERNAL_OBJECTS =

src/IO/FILE_READERS/test_logreader: src/IO/FILE_READERS/CMakeFiles/test_logreader.dir/test_logreader.cpp.o
src/IO/FILE_READERS/test_logreader: src/IO/FILE_READERS/CMakeFiles/test_logreader.dir/build.make
src/IO/FILE_READERS/test_logreader: src/IO/FILE_READERS/libFILE_READERS_SRC.a
src/IO/FILE_READERS/test_logreader: src/MATRIX/libMATRIX_SRC.a
src/IO/FILE_READERS/test_logreader: src/STRING_SUPPORT/libSTRING_SUPPORT_SRC.a
src/IO/FILE_READERS/test_logreader: src/IO/FILE_READERS/CMakeFiles/test_logreader.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX executable test_logreader"
	cd /home/x_nicro/LOE-CTP-FRAG/QC_Tools/src/IO/FILE_READERS && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/test_logreader.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/IO/FILE_READERS/CMakeFiles/test_logreader.dir/build: src/IO/FILE_READERS/test_logreader
.PHONY : src/IO/FILE_READERS/CMakeFiles/test_logreader.dir/build

src/IO/FILE_READERS/CMakeFiles/test_logreader.dir/requires: src/IO/FILE_READERS/CMakeFiles/test_logreader.dir/test_logreader.cpp.o.requires
.PHONY : src/IO/FILE_READERS/CMakeFiles/test_logreader.dir/requires

src/IO/FILE_READERS/CMakeFiles/test_logreader.dir/clean:
	cd /home/x_nicro/LOE-CTP-FRAG/QC_Tools/src/IO/FILE_READERS && $(CMAKE_COMMAND) -P CMakeFiles/test_logreader.dir/cmake_clean.cmake
.PHONY : src/IO/FILE_READERS/CMakeFiles/test_logreader.dir/clean

src/IO/FILE_READERS/CMakeFiles/test_logreader.dir/depend:
	cd /home/x_nicro/LOE-CTP-FRAG/QC_Tools && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/x_nicro/LOE-CTP-FRAG/QC_Tools /home/x_nicro/LOE-CTP-FRAG/QC_Tools/src/IO/FILE_READERS /home/x_nicro/LOE-CTP-FRAG/QC_Tools /home/x_nicro/LOE-CTP-FRAG/QC_Tools/src/IO/FILE_READERS /home/x_nicro/LOE-CTP-FRAG/QC_Tools/src/IO/FILE_READERS/CMakeFiles/test_logreader.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/IO/FILE_READERS/CMakeFiles/test_logreader.dir/depend
