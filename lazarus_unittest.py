##############################################################################
#
# Lazarus Unit Tests
#
# This test file is deprecated!
#
# Author: Victor Hanson-Smith
# Contact: victorhs@cs.uoregon.edu
#
###############################################################################
from includeMe import *
from utest_tree import *
from utest_engine import *

try:
    test_tree()
    t = testEngine()
    t.test_engine()
except AssertionError:
    print "\n\nThe unit tests failed.  Please review the previous error messages for more information."
    exit(1)
print "\n\nSUCCESS! All tests returned successful results.\n"