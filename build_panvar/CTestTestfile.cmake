# CMake generated Testfile for 
# Source directory: /Users/davide.bolognini/panvar
# Build directory: /Users/davide.bolognini/panvar/build_panvar
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(panvar_smoke_c4 "/usr/bin/env" "bash" "/Users/davide.bolognini/panvar/tests/smoke.sh" "/Users/davide.bolognini/panvar/build_panvar/panvar" "/Users/davide.bolognini/panvar/tests/real_data/c4.gfa" "/Users/davide.bolognini/panvar/build_panvar/smoke_c4")
set_tests_properties(panvar_smoke_c4 PROPERTIES  _BACKTRACE_TRIPLES "/Users/davide.bolognini/panvar/CMakeLists.txt;97;add_test;/Users/davide.bolognini/panvar/CMakeLists.txt;0;")
