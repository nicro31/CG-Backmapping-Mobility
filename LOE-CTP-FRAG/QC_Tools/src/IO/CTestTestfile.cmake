# CMake generated Testfile for 
# Source directory: /home/x_nicro/LOE-CTP-FRAG/QC_Tools/src/IO
# Build directory: /home/x_nicro/LOE-CTP-FRAG/QC_Tools/src/IO
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
ADD_TEST(test_io "test_io" "-p_P" "../../../GAUSSIANFILES/30/30_pair.pun" "-p_1" "../../../GAUSSIANFILES/30/ref.pun" "-p_2" "../../../GAUSSIANFILES/30/30_2.pun")
ADD_TEST(test_argumentparser "test_argumentparser")
ADD_TEST(test_io_script "/bin/bash" "/home/x_nicro/LOE-CTP-FRAG/QC_Tools/src/IO/test_script_io.sh" "/home/x_nicro/LOE-CTP-FRAG/QC_Tools/src/IO")
SUBDIRS(FILE_READERS)
SUBDIRS(ARGUMENTS)
