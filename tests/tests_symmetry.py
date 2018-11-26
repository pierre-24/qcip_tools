from tests import QcipToolsTestCase
from qcip_tools import symmetry


class SymmetryTestCase(QcipToolsTestCase):

    def test_symmetry_element(self):
        s = symmetry.SymmetryElement.S(4, 7)
        print(s)
