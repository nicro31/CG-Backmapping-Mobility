# CMake generated Testfile for 
# Source directory: /home/x_nicro/LOE-CTP-FRAG/Approx_MVBB/ApproxMVBB/tests
# Build directory: /home/x_nicro/LOE-CTP-FRAG/Approx_MVBB/ApproxMVBB/ApproxMVBB_Build/tests
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(ApproxMVBBTest-ConvexHull "/home/x_nicro/LOE-CTP-FRAG/Approx_MVBB/ApproxMVBB/ApproxMVBB_Build/tests/bin/ApproxMVBBTest-ConvexHull")
add_test(ApproxMVBBTest-MinAreaRect "/home/x_nicro/LOE-CTP-FRAG/Approx_MVBB/ApproxMVBB/ApproxMVBB_Build/tests/bin/ApproxMVBBTest-MinAreaRect")
add_test(ApproxMVBBTest-Diameter "/home/x_nicro/LOE-CTP-FRAG/Approx_MVBB/ApproxMVBB/ApproxMVBB_Build/tests/bin/ApproxMVBBTest-Diameter")
add_test(ApproxMVBBTest-DiameterOOBB "/home/x_nicro/LOE-CTP-FRAG/Approx_MVBB/ApproxMVBB/ApproxMVBB_Build/tests/bin/ApproxMVBBTest-DiameterOOBB")
add_test(ApproxMVBBTest-MVBB "/home/x_nicro/LOE-CTP-FRAG/Approx_MVBB/ApproxMVBB/ApproxMVBB_Build/tests/bin/ApproxMVBBTest-MVBB")
subdirs("../../../../thirdparty/googletest-build")
