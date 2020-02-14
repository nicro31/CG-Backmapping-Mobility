# CMake generated Testfile for 
# Source directory: /home/x_nicro/LOE-CTP-FRAG/QC_Tools/src/IO
# Build directory: /home/x_nicro/LOE-CTP-FRAG/QC_Tools/QC_Tools_Build/src/IO
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(test_io "test_io" "-p_P" "../../../GAUSSIANFILES/30/30_pair.pun" "-p_1" "../../../GAUSSIANFILES/30/ref.pun" "-p_2" "../../../GAUSSIANFILES/30/30_2.pun")
add_test(test_argumentparser "test_argumentparser")
add_test(test_io_script "/bin/bash" "/home/x_nicro/LOE-CTP-FRAG/QC_Tools/src/IO/test_script_io.sh" "/home/x_nicro/LOE-CTP-FRAG/QC_Tools/src/IO")
subdirs("FILE_READERS")
subdirs("ARGUMENTS")
