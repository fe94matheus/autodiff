# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.27

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/local/bin/cmake

# The command to remove a file.
RM = /usr/local/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/felipe/Downloads/autodiff

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/felipe/Downloads/autodiff/build

# Include any dependencies generated for this target.
include CMakeFiles/AutoDiff.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/AutoDiff.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/AutoDiff.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/AutoDiff.dir/flags.make

CMakeFiles/AutoDiff.dir/src/autodiff.cpp.o: CMakeFiles/AutoDiff.dir/flags.make
CMakeFiles/AutoDiff.dir/src/autodiff.cpp.o: /home/felipe/Downloads/autodiff/src/autodiff.cpp
CMakeFiles/AutoDiff.dir/src/autodiff.cpp.o: CMakeFiles/AutoDiff.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/felipe/Downloads/autodiff/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/AutoDiff.dir/src/autodiff.cpp.o"
	/home/felipe/gcc-13-build/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/AutoDiff.dir/src/autodiff.cpp.o -MF CMakeFiles/AutoDiff.dir/src/autodiff.cpp.o.d -o CMakeFiles/AutoDiff.dir/src/autodiff.cpp.o -c /home/felipe/Downloads/autodiff/src/autodiff.cpp

CMakeFiles/AutoDiff.dir/src/autodiff.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/AutoDiff.dir/src/autodiff.cpp.i"
	/home/felipe/gcc-13-build/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/felipe/Downloads/autodiff/src/autodiff.cpp > CMakeFiles/AutoDiff.dir/src/autodiff.cpp.i

CMakeFiles/AutoDiff.dir/src/autodiff.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/AutoDiff.dir/src/autodiff.cpp.s"
	/home/felipe/gcc-13-build/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/felipe/Downloads/autodiff/src/autodiff.cpp -o CMakeFiles/AutoDiff.dir/src/autodiff.cpp.s

CMakeFiles/AutoDiff.dir/src/main.cpp.o: CMakeFiles/AutoDiff.dir/flags.make
CMakeFiles/AutoDiff.dir/src/main.cpp.o: /home/felipe/Downloads/autodiff/src/main.cpp
CMakeFiles/AutoDiff.dir/src/main.cpp.o: CMakeFiles/AutoDiff.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/felipe/Downloads/autodiff/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/AutoDiff.dir/src/main.cpp.o"
	/home/felipe/gcc-13-build/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/AutoDiff.dir/src/main.cpp.o -MF CMakeFiles/AutoDiff.dir/src/main.cpp.o.d -o CMakeFiles/AutoDiff.dir/src/main.cpp.o -c /home/felipe/Downloads/autodiff/src/main.cpp

CMakeFiles/AutoDiff.dir/src/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/AutoDiff.dir/src/main.cpp.i"
	/home/felipe/gcc-13-build/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/felipe/Downloads/autodiff/src/main.cpp > CMakeFiles/AutoDiff.dir/src/main.cpp.i

CMakeFiles/AutoDiff.dir/src/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/AutoDiff.dir/src/main.cpp.s"
	/home/felipe/gcc-13-build/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/felipe/Downloads/autodiff/src/main.cpp -o CMakeFiles/AutoDiff.dir/src/main.cpp.s

# Object files for target AutoDiff
AutoDiff_OBJECTS = \
"CMakeFiles/AutoDiff.dir/src/autodiff.cpp.o" \
"CMakeFiles/AutoDiff.dir/src/main.cpp.o"

# External object files for target AutoDiff
AutoDiff_EXTERNAL_OBJECTS =

AutoDiff: CMakeFiles/AutoDiff.dir/src/autodiff.cpp.o
AutoDiff: CMakeFiles/AutoDiff.dir/src/main.cpp.o
AutoDiff: CMakeFiles/AutoDiff.dir/build.make
AutoDiff: CMakeFiles/AutoDiff.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --bold --progress-dir=/home/felipe/Downloads/autodiff/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Linking CXX executable AutoDiff"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/AutoDiff.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/AutoDiff.dir/build: AutoDiff
.PHONY : CMakeFiles/AutoDiff.dir/build

CMakeFiles/AutoDiff.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/AutoDiff.dir/cmake_clean.cmake
.PHONY : CMakeFiles/AutoDiff.dir/clean

CMakeFiles/AutoDiff.dir/depend:
	cd /home/felipe/Downloads/autodiff/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/felipe/Downloads/autodiff /home/felipe/Downloads/autodiff /home/felipe/Downloads/autodiff/build /home/felipe/Downloads/autodiff/build /home/felipe/Downloads/autodiff/build/CMakeFiles/AutoDiff.dir/DependInfo.cmake "--color=$(COLOR)"
.PHONY : CMakeFiles/AutoDiff.dir/depend

