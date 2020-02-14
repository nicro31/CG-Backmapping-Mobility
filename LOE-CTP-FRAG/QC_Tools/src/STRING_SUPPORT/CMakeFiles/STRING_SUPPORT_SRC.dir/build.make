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
include src/STRING_SUPPORT/CMakeFiles/STRING_SUPPORT_SRC.dir/depend.make

# Include the progress variables for this target.
include src/STRING_SUPPORT/CMakeFiles/STRING_SUPPORT_SRC.dir/progress.make

# Include the compile flags for this target's objects.
include src/STRING_SUPPORT/CMakeFiles/STRING_SUPPORT_SRC.dir/flags.make

src/STRING_SUPPORT/CMakeFiles/STRING_SUPPORT_SRC.dir/string_support.cpp.o: src/STRING_SUPPORT/CMakeFiles/STRING_SUPPORT_SRC.dir/flags.make
src/STRING_SUPPORT/CMakeFiles/STRING_SUPPORT_SRC.dir/string_support.cpp.o: src/STRING_SUPPORT/string_support.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/x_nicro/LOE-CTP-FRAG/QC_Tools/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object src/STRING_SUPPORT/CMakeFiles/STRING_SUPPORT_SRC.dir/string_support.cpp.o"
	cd /home/x_nicro/LOE-CTP-FRAG/QC_Tools/src/STRING_SUPPORT && /software/sse/wrappers/icpc   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/STRING_SUPPORT_SRC.dir/string_support.cpp.o -c /home/x_nicro/LOE-CTP-FRAG/QC_Tools/src/STRING_SUPPORT/string_support.cpp

src/STRING_SUPPORT/CMakeFiles/STRING_SUPPORT_SRC.dir/string_support.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/STRING_SUPPORT_SRC.dir/string_support.cpp.i"
	cd /home/x_nicro/LOE-CTP-FRAG/QC_Tools/src/STRING_SUPPORT && /software/sse/wrappers/icpc  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/x_nicro/LOE-CTP-FRAG/QC_Tools/src/STRING_SUPPORT/string_support.cpp > CMakeFiles/STRING_SUPPORT_SRC.dir/string_support.cpp.i

src/STRING_SUPPORT/CMakeFiles/STRING_SUPPORT_SRC.dir/string_support.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/STRING_SUPPORT_SRC.dir/string_support.cpp.s"
	cd /home/x_nicro/LOE-CTP-FRAG/QC_Tools/src/STRING_SUPPORT && /software/sse/wrappers/icpc  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/x_nicro/LOE-CTP-FRAG/QC_Tools/src/STRING_SUPPORT/string_support.cpp -o CMakeFiles/STRING_SUPPORT_SRC.dir/string_support.cpp.s

src/STRING_SUPPORT/CMakeFiles/STRING_SUPPORT_SRC.dir/string_support.cpp.o.requires:
.PHONY : src/STRING_SUPPORT/CMakeFiles/STRING_SUPPORT_SRC.dir/string_support.cpp.o.requires

src/STRING_SUPPORT/CMakeFiles/STRING_SUPPORT_SRC.dir/string_support.cpp.o.provides: src/STRING_SUPPORT/CMakeFiles/STRING_SUPPORT_SRC.dir/string_support.cpp.o.requires
	$(MAKE) -f src/STRING_SUPPORT/CMakeFiles/STRING_SUPPORT_SRC.dir/build.make src/STRING_SUPPORT/CMakeFiles/STRING_SUPPORT_SRC.dir/string_support.cpp.o.provides.build
.PHONY : src/STRING_SUPPORT/CMakeFiles/STRING_SUPPORT_SRC.dir/string_support.cpp.o.provides

src/STRING_SUPPORT/CMakeFiles/STRING_SUPPORT_SRC.dir/string_support.cpp.o.provides.build: src/STRING_SUPPORT/CMakeFiles/STRING_SUPPORT_SRC.dir/string_support.cpp.o

# Object files for target STRING_SUPPORT_SRC
STRING_SUPPORT_SRC_OBJECTS = \
"CMakeFiles/STRING_SUPPORT_SRC.dir/string_support.cpp.o"

# External object files for target STRING_SUPPORT_SRC
STRING_SUPPORT_SRC_EXTERNAL_OBJECTS =

src/STRING_SUPPORT/libSTRING_SUPPORT_SRC.a: src/STRING_SUPPORT/CMakeFiles/STRING_SUPPORT_SRC.dir/string_support.cpp.o
src/STRING_SUPPORT/libSTRING_SUPPORT_SRC.a: src/STRING_SUPPORT/CMakeFiles/STRING_SUPPORT_SRC.dir/build.make
src/STRING_SUPPORT/libSTRING_SUPPORT_SRC.a: src/STRING_SUPPORT/CMakeFiles/STRING_SUPPORT_SRC.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX static library libSTRING_SUPPORT_SRC.a"
	cd /home/x_nicro/LOE-CTP-FRAG/QC_Tools/src/STRING_SUPPORT && $(CMAKE_COMMAND) -P CMakeFiles/STRING_SUPPORT_SRC.dir/cmake_clean_target.cmake
	cd /home/x_nicro/LOE-CTP-FRAG/QC_Tools/src/STRING_SUPPORT && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/STRING_SUPPORT_SRC.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/STRING_SUPPORT/CMakeFiles/STRING_SUPPORT_SRC.dir/build: src/STRING_SUPPORT/libSTRING_SUPPORT_SRC.a
.PHONY : src/STRING_SUPPORT/CMakeFiles/STRING_SUPPORT_SRC.dir/build

src/STRING_SUPPORT/CMakeFiles/STRING_SUPPORT_SRC.dir/requires: src/STRING_SUPPORT/CMakeFiles/STRING_SUPPORT_SRC.dir/string_support.cpp.o.requires
.PHONY : src/STRING_SUPPORT/CMakeFiles/STRING_SUPPORT_SRC.dir/requires

src/STRING_SUPPORT/CMakeFiles/STRING_SUPPORT_SRC.dir/clean:
	cd /home/x_nicro/LOE-CTP-FRAG/QC_Tools/src/STRING_SUPPORT && $(CMAKE_COMMAND) -P CMakeFiles/STRING_SUPPORT_SRC.dir/cmake_clean.cmake
.PHONY : src/STRING_SUPPORT/CMakeFiles/STRING_SUPPORT_SRC.dir/clean

src/STRING_SUPPORT/CMakeFiles/STRING_SUPPORT_SRC.dir/depend:
	cd /home/x_nicro/LOE-CTP-FRAG/QC_Tools && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/x_nicro/LOE-CTP-FRAG/QC_Tools /home/x_nicro/LOE-CTP-FRAG/QC_Tools/src/STRING_SUPPORT /home/x_nicro/LOE-CTP-FRAG/QC_Tools /home/x_nicro/LOE-CTP-FRAG/QC_Tools/src/STRING_SUPPORT /home/x_nicro/LOE-CTP-FRAG/QC_Tools/src/STRING_SUPPORT/CMakeFiles/STRING_SUPPORT_SRC.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/STRING_SUPPORT/CMakeFiles/STRING_SUPPORT_SRC.dir/depend

