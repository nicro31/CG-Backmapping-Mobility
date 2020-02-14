# CMake generated Testfile for 
# Source directory: /home/x_nicro/LOE-CTP-FRAG/QC_Tools
# Build directory: /home/x_nicro/LOE-CTP-FRAG/QC_Tools/QC_Tools_Build
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(test_calc_J_script "/bin/bash" "/home/x_nicro/LOE-CTP-FRAG/QC_Tools/src/test_script_calc_J.sh" "/home/x_nicro/LOE-CTP-FRAG")
subdirs("src/MATRIX")
subdirs("src/STRING_SUPPORT")
subdirs("src/IO")
subdirs("src/PARAMETERS")
subdirs("src/QC_FUNCTIONS")
