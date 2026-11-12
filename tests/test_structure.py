import signal
import tempfile
import unittest
from contextlib import contextmanager
from pathlib import Path

from drawtetrado import structure


ROOT = Path(__file__).resolve().parents[1]
DATA = Path(__file__).resolve().parent / "data"


@contextmanager
def time_limit(seconds):
    """Fail instead of hanging if a regression brings back an infinite loop."""
    if not hasattr(signal, "SIGALRM"):
        yield
        return

    def handler(signum, frame):
        raise TimeoutError("test exceeded {} seconds".format(seconds))

    previous = signal.signal(signal.SIGALRM, handler)
    signal.setitimer(signal.ITIMER_REAL, seconds)
    try:
        yield
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)
        signal.signal(signal.SIGALRM, previous)


def complement(lw):
    return "cHW" if lw == "cWH" else "cWH"


def flip_sign(sign):
    return "-" if sign == "+" else "+"


def pair_lw(structure_obj, nt1, nt2):
    """Return the Leontis-Westhof label of nt1 -> nt2."""
    for pair in structure_obj.basePairs:
        lw = pair.get("lw")
        if lw not in ("cWH", "cHW"):
            continue
        if pair["nt1"] == nt1 and pair["nt2"] == nt2:
            return lw
        if pair["nt1"] == nt2 and pair["nt2"] == nt1:
            return complement(lw)
    return None


def tetrad_onz(structure_obj, quadruplex_id, names):
    names = set(names)
    for tetrad in structure_obj.tetrads[quadruplex_id].values():
        bases = {tetrad["nt1"], tetrad["nt2"], tetrad["nt3"], tetrad["nt4"]}
        if names <= bases:
            return tetrad["onz"]
    raise AssertionError("no tetrad contains {}".format(sorted(names)))


def position_pairing_type(quadruplex, structure_obj, tetrad):
    """Visual cWX label of the tetrad edges in increasing position order."""
    onz = tetrad_onz(structure_obj, 0, tetrad)
    pairing = "cWH" if onz.endswith("+") else "cHW"
    if quadruplex.nucl_quad[tetrad[0]].mirrored:
        pairing = complement(pairing)
    return pairing


def visual_pair_lw(quadruplex, structure_obj, tetrad, nt1, nt2):
    """Visual label of nt1 -> nt2 as drawn in the tetrad."""
    index_1 = tetrad.index(nt1)
    index_2 = tetrad.index(nt2)
    pairing = position_pairing_type(quadruplex, structure_obj, tetrad)

    if (index_2 - index_1) % 4 == 1:
        return pairing
    if (index_1 - index_2) % 4 == 1:
        return complement(pairing)
    raise AssertionError("{} and {} are not adjacent".format(nt1, nt2))


class TetradOrientationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.structure = structure.Structure().fromFile(str(DATA / "8taa.json"))

    def test_prepared_order_follows_lw(self):
        quadruplex = structure.Quadruplex(self.structure, 0)
        for tetrad in quadruplex.tetrads:
            onz = tetrad_onz(self.structure, 0, tetrad)
            expected = "cWH" if onz.endswith("+") else "cHW"
            for index in range(4):
                nt1 = tetrad[index]
                nt2 = tetrad[(index + 1) % 4]
                self.assertEqual(pair_lw(self.structure, nt1, nt2), expected)

    def test_same_rotations_share_alignment_group(self):
        quadruplex = structure.Quadruplex(self.structure, 0)
        self.assertEqual(quadruplex.GetSameRotations(), [0, 0, 0, 0])

    def test_reflected_tetrads_use_mirrored_geometry(self):
        quadruplex = structure.Quadruplex(self.structure, 0)
        quadruplex.Optimize()

        mirrored = {
            tetrad[0]: quadruplex.nucl_quad[tetrad[0]].mirrored
            for tetrad in quadruplex.tetrads
        }
        self.assertTrue(any(mirrored.values()))

        for tetrad in quadruplex.tetrads:
            for name in tetrad:
                nucl = quadruplex.nucl_quad[name]
                expected = nucl.GetOnzPlusMinus()
                if nucl.mirrored:
                    expected = flip_sign(expected)
                self.assertEqual(nucl.GetGeometrySign(), expected)

    def test_optimizer_preserves_lw_pairing(self):
        quadruplex = structure.Quadruplex(self.structure, 0)
        quadruplex.Optimize()

        for tetrad in quadruplex.tetrads:
            pairing = position_pairing_type(quadruplex, self.structure, tetrad)
            for index in range(4):
                nt1 = tetrad[index]
                nt2 = tetrad[(index + 1) % 4]
                self.assertEqual(pair_lw(self.structure, nt1, nt2), pairing)

        tetrad = next(
            tetrad
            for tetrad in quadruplex.tetrads
            if "A.DG21" in tetrad and "A.DG24" in tetrad
        )
        self.assertEqual(
            visual_pair_lw(quadruplex, self.structure, tetrad, "A.DG21", "A.DG24"),
            "cWH",
        )


class TractAlignmentTests(unittest.TestCase):
    def assert_continuous_tracts(self, fixture):
        structure_obj = structure.Structure().fromFile(str(DATA / fixture))
        quadruplex = structure.Quadruplex(structure_obj, 0)
        quadruplex.Optimize()
        tract_of = quadruplex.GetTractOf()

        for tract_index in range(4):
            positions = {}
            for name, (group, index) in tract_of.items():
                if group == 0 and index == tract_index:
                    positions[name] = quadruplex.nucl_quad[name].position

            self.assertGreater(len(positions), 1)
            self.assertEqual(
                len(set(positions.values())),
                1,
                msg="tract {} is not one column: {}".format(tract_index, positions),
            )

        return structure_obj, quadruplex

    def test_7ys7_tracts_are_columns(self):
        structure_obj, quadruplex = self.assert_continuous_tracts("7ys7.json")

        tract = structure_obj.tracts[0][0][0]
        self.assertEqual(tract, ["A.DG2", "A.DG1", "A.DG16"])
        column = quadruplex.nucl_quad[tract[0]].position
        for name in tract:
            self.assertEqual(quadruplex.nucl_quad[name].position, column)
        self.assertEqual(
            [quadruplex.nucl_quad[name].tetrade_no for name in tract],
            [2, 1, 0],
        )

    def test_8taa_tracts_are_columns(self):
        self.assert_continuous_tracts("8taa.json")


class LegacyInputTests(unittest.TestCase):
    def test_legacy_input_without_lw_does_not_hang(self):
        structure_obj = structure.Structure().fromFile(str(ROOT / "examples" / "1myq.json"))
        with time_limit(10):
            quadruplex = structure.Quadruplex(structure_obj, 0)

        self.assertGreater(len(quadruplex.tetrads), 0)
        for tetrad in quadruplex.tetrads:
            self.assertEqual(len(tetrad), 4)


class DrawingTests(unittest.TestCase):
    def test_draw_8taa(self):
        from drawtetrado.main import DrawFromFile
        from drawtetrado.svg_painter import Config

        with tempfile.TemporaryDirectory() as tmp:
            prefix = str(Path(tmp) / "8taa")
            DrawFromFile(str(DATA / "8taa.json"), prefix, Config(1.0))

            svg_path = Path(prefix + "_0.svg")
            self.assertTrue(svg_path.exists())

            svg = svg_path.read_text()
            self.assertIn("A.DG21", svg)
            self.assertIn("A.DG24", svg)


if __name__ == "__main__":
    unittest.main()
