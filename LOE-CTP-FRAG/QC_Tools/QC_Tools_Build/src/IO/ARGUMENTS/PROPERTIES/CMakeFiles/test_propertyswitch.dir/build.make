# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.12

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


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
CMAKE_COMMAND = /software/sse/manual/CMake/3.12.1/bin/cmake

# The command to remove a file.
RM = /software/sse/manual/CMake/3.12.1/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/x_nicro/LOE-CTP-FRAG

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/x_nicro/LOE-CTP-FRAG

# Include any dependencies generated for this target.
include QC_Tools/QC_Tools_Build/src/IO/ARGUMENTS/PROPERTIES/CMakeFiles/test_propertyswitch.dir/depend.make

# Include the progress variables for this target.
include QC_Tools/QC_Tools_Build/src/IO/ARGUMENTS/PROPERTIES/CMakeFiles/test_propertyswitch.dir/progress.make

# Include the compile flags for this target's objects.
include QC_Tools/QC_Tools_Build/src/IO/ARGUMENTS/PROPERTIES/CMakeFiles/test_propertyswitch.dir/flags.make

QC_Tools/QC_Tools_Build/src/IO/ARGUMENTS/PROPERTIES/CMakeFiles/test_propertyswitch.dir/test_propertyswitch.cpp.o: QC_Tools/QC_Tools_Build/src/IO/ARGUMENTS/PROPERTIES/CMakeFiles/test_propertyswitch.dir/flags.make
QC_Tools/QC_Tools_Build/src/IO/ARGUMENTS/PROPERTIES/CMakeFiles/test_propertyswitch.dir/test_propertyswitch.cpp.o: QC_Tools/src/IO/ARGUMENTS/PROPERTIES/test_propertyswitch.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/x_nicro/LOE-CTP-FRAG/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object QC_Tools/QC_Tools_Build/src/IO/ARGUMENTS/PROPERTIES/CMakeFiles/test_propertyswitch.dir/test_propertyswitch.cpp.o"
	cd /home/x_nicro/LOE-CTP-FRAG/QC_Tools/QC_Tools_Build/src/IO/ARGUMENTS/PROPERTIES && /software/sse/manual/mpprun/4.1.2/nsc-wrappers/icpc  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/test_propertyswitch.dir/test_propertyswitch.cpp.o -c /home/x_nicro/LOE-CTP-FRAG/QC_Tools/src/IO/ARGUMENTS/PROPERTIES/test_propertyswitch.cpp

QC_Tools/QC_Tools_Build/src/IO/ARGUMENTS/PROPERTIES/CMakeFiles/test_propertyswitch.dir/test_propertyswitch.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/test_propertyswitch.dir/test_propertyswitch.cpp.i"
	cd /home/x_nicro/LOE-CTP-FRAG/QC_Tools/QC_Tools_Build/src/IO/ARGUMENTS/PROPERTIES && /software/sse/manual/mpprun/4.1.2/nsc-wrappers/icpc $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/x_nicro/LOE-CTP-FRAG/QC_Tools/src/IO/ARGUMENTS/PROPERTIES/test_propertyswitch.cpp > CMakeFiles/test_propertyswitch.dir/test_propertyswitch.cpp.i

QC_Tools/QC_Tools_Build/src/IO/ARGUMENTS/PROPERTIES/CMakeFiles/test_propertyswitch.dir/test_propertyswitch.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/test_propertyswitch.dir/test_propertyswitch.cpp.s"
	cd /home/x_nicro/LOE-CTP-FRAG/QC_Tools/QC_Tools_Build/src/IO/ARGUMENTS/PROPERTIES && /software/sse/manual/mpprun/4.1.2/nsc-wrappers/icpc $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/x_nicro/LOE-CTP-FRAG/QC_Tools/src/IO/ARGUMENTS/PROPERTIES/test_propertyswitch.cpp -o CMakeFiles/test_propertyswitch.dir/test_propertyswitch.cpp.s

# Object files for target test_propertyswitch
test_propertyswitch_OBJECTS = \
"CMakeFiles/test_propertyswitch.dir/test_propertyswitch.cpp.o"

# External object files for target test_propertyswitch
test_propertyswitch_EXTERNAL_OBJECTS =

QC_Tools/QC_Tools_Build/src/IO/ARGUMENTS/PROPERTIES/bin/test_propertyswitch: QC_Tools/QC_Tools_Build/src/IO/ARGUMENTS/PROPERTIES/CMakeFiles/test_propertyswitch.dir/test_propertyswitch.cpp.o
QC_Tools/QC_Tools_Build/src/IO/ARGUMENTS/PROPERTIES/bin/test_propertyswitch: QC_Tools/QC_Tools_Build/src/IO/ARGUMENTS/PROPERTIES/CMakeFiles/test_propertyswitch.dir/build.make
QC_Tools/QC_Tools_Build/src/IO/ARGUMENTS/PROPERTIES/bin/test_propertyswitch: QC_Tools/QC_Tools_Build/src/IO/ARGUMENTS/PROPERTIES/libPROPERTIES_SRC.a
QC_Tools/QC_Tools_Build/src/IO/ARGUMENTS/PROPERTIES/bin/test_propertyswitch: QC_Tools/QC_Tools_Build/src/STRING_SUPPORT/libSTRING_SUPPORT_SRC.a
QC_Tools/QC_Tools_Build/src/IO/ARGUMENTS/PROPERTIES/bin/test_propertyswitch: QC_Tools/QC_Tools_Build/src/IO/ARGUMENTS/PROPERTIES/CMakeFiles/test_propertyswitch.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/x_nicro/LOE-CTP-FRAG/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable bin/test_propertyswitch"
	cd /home/x_nicro/LOE-CTP-FRAG/QC_Tools/QC_Tools_Build/src/IO/ARGUMENTS/PROPERTIES && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/test_propertyswitch.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
QC_Tools/QC_Tools_Build/src/IO/ARGUMENTS/PROPERTIES/CMakeFiles/test_propertyswitch.dir/build: QC_Tools/QC_Tools_Build/src/IO/ARGUMENTS/PROPERTIES/bin/test_propertyswitch

.PHONY : QC_Tools/QC_Tools_Build/src/IO/ARGUMENTS/PROPERTIES/CMakeFiles/test_propertyswitch.dir/build

QC_Tools/QC_Tools_Build/src/IO/ARGUMENTS/PROPERTIES/CMakeFiles/test_propertyswitch.dir/clean:
	cd /home/x_nicro/LOE-CTP-FRAG/QC_Tools/QC_Tools_Build/src/IO/ARGUMENTS/PROPERTIES && $(CMAKE_COMMAND) -P CMakeFiles/test_propertyswitch.dir/cmake_clean.cmake
.PHONY : QC_Tools/QC_Tools_Build/src/IO/ARGUMENTS/PROPERTIES/CMakeFiles/test_propertyswitch.dir/clean

QC_Tools/QC_Tools_Build/src/IO/ARGUMENTS/PROPERTIES/CMakeFiles/test_propertyswitch.dir/depend:
	cd /home/x_nicro/LOE-CTP-FRAG && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/x_nicro/LOE-CTP-FRAG /home/x_nicro/LOE-CTP-FRAG/QC_Tools/src/IO/ARGUMENTS/PROPERTIES /home/x_nicro/LOE-CTP-FRAG /home/x_nicro/LOE-CTP-FRAG/QC_Tools/QC_Tools_Build/src/IO/ARGUMENTS/PROPERTIES /home/x_nicro/LOE-CTP-FRAG/QC_Tools/QC_Tools_Build/src/IO/ARGUMENTS/PROPERTIES/CMakeFiles/test_propertyswitch.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : QC_Tools/QC_Tools_Build/src/IO/ARGUMENTS/PROPERTIES/CMakeFiles/test_propertyswitch.dir/depend

