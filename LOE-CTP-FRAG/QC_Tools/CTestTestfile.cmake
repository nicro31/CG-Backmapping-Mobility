# CMake generated Testfile for 
# Source directory: /home/x_nicro/LOE-CTP-FRAG/QC_Tools
# Build directory: /home/x_nicro/LOE-CTP-FRAG/QC_Tools
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
ADD_TEST(test_calc_J_script "/bin/bash" "/home/x_nicro/LOE-CTP-FRAG/QC_Tools/src/test_script_calc_J.sh" "/home/x_nicro/LOE-CTP-FRAG/QC_Tools")
SUBDIRS(src/MATRIX)
SUBDIRS(src/STRING_SUPPORT)
SUBDIRS(src/IO)
SUBDIRS(src/PARAMETERS)
SUBDIRS(src/QC_FUNCTIONS)
