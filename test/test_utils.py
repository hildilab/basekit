import unittest

from basekit.utils import (
    try_int
)


class UtilsTestCase( unittest.TestCase ):
    def test_try_int( self ):
        self.assertEqual( 3, try_int("3") )
        


