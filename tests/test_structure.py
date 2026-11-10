import signal
import tempfile
import unittest
from contextlib import contextmanager
from pathlib import Path

from drawtetrado import structure


ROOT = Path(__file__).resolve().parents[1]
DATA = Path(__file__).resolve().parent / "data"
ROTATIONS = {
    (0, 1, 2, 3),
    (3, 0, 1, 2),
    (2, 3, 0, 1),
    (1, 2, 3, 0),
}


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


def pair_lw(structure_obj, nt1, nt2):
    """Return the Leontis-Westhof label of nt1 -> nt2."""
    for pair in structure_obj.basePairs:
        lw = pair.get("lw")
        if lw not in ("cWH", "cHW"):
            continue
        if pair["nt1"] == nt1 and pair["nt2"] == nt2:
            return lw
        if pair["nt1"] == nt2 and pair["nt2"] == nt1:
            return "cHW" if lw == "cWH" else "cWH"
    return None


def tetrad_onz(structure_obj, quadruplex_id, names):
    names = set(names)
    for tetrad in structure_obj.tetrads[quadruplex_id].values():
        bases = {tetrad["nt1"], tetrad["nt2"], tetrad["nt3"], tetrad["nt4"]}
        if names <= bases:
            return tetrad["onz"]
    raise AssertionError("no tetrad contains {}".format(sorted(names)))


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

    def test_same_rotations_separate_polarity_blocks(self):
        quadruplex = structure.Quadruplex(self.structure, 0)
        self.assertEqual(quadruplex.GetSameRotations(), [0, 0, 1, 1])

    def test_optimizer_uses_rotations_only(self):
        quadruplex = structure.Quadruplex(self.structure, 0)
        prepared = [list(tetrad) for tetrad in quadruplex.tetrads]

        quadruplex.Optimize()

        for before, after in zip(prepared, quadruplex.tetrads):
            rotations = {
                tuple(before[index:] + before[:index]) for index in range(4)
            }
            self.assertIn(tuple(after), rotations, msg=str(after))

    def test_optimizer_never_searches_reflections(self):
        import optimizer

        quadruplex = structure.Quadruplex(self.structure, 0)
        edges = quadruplex.GetNucleotidesPositions()
        alignments = quadruplex.GetAlignments()
        # Force all levels to share one alignment group. This is only
        # satisfiable with reflections, so the rotation-only optimizer must
        # fall back to free rotations.
        rotations = [0] * len(quadruplex.tetrads)

        optimized = optimizer.solve(edges, rotations, alignments)

        for level in range(len(quadruplex.tetrads)):
            positions = tuple(optimized[level * 4:(level + 1) * 4])
            self.assertIn(positions, ROTATIONS)

    def test_optimizer_preserves_lw_pairing(self):
        quadruplex = structure.Quadruplex(self.structure, 0)
        quadruplex.Optimize()

        for tetrad in quadruplex.tetrads:
            onz = tetrad_onz(self.structure, 0, tetrad)
            expected = "cWH" if onz.endswith("+") else "cHW"
            for index in range(4):
                nt1 = tetrad[index]
                nt2 = tetrad[(index + 1) % 4]
                self.assertEqual(pair_lw(self.structure, nt1, nt2), expected)

        tetrad = next(
            tetrad
            for tetrad in quadruplex.tetrads
            if "A.DG21" in tetrad and "A.DG24" in tetrad
        )
        index_21 = tetrad.index("A.DG21")
        index_24 = tetrad.index("A.DG24")
        self.assertEqual((index_24 - index_21) % 4, 1)
        self.assertEqual(pair_lw(self.structure, "A.DG21", "A.DG24"), "cWH")


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
