import unittest
from diagram import generate_frequencies


class TestGenerateDiagrams(unittest.TestCase):
    def test_generate_frequencies(self):
        k1 = [1, 1]
        freqs_k1 = [[1, 0, 0], [-1, 0, 0]]
        self.assertEqual(generate_frequencies(k1), freqs_k1)

        k2 = [0, 1]
        freqs_k2 = [[1, 1, 0], [1, -1, 0], [-1, 1, 0], [-1, -1, 0]]
        self.assertEqual(generate_frequencies(k2), freqs_k2)

        k2p = [1, 0]
        freqs_k2p = [[1, 0, 1], [1, 0, -1], [-1, 0, 1], [-1, 0, -1]]
        self.assertEqual(generate_frequencies(k2p), freqs_k2p)

        k3 = [0, 0]
        freqs_k3 = [[1, 1, 1], [1, 1, -1], [1, -1, 1], [1, -1, -1],
                    [-1, 1, 1], [-1, 1, -1], [-1, -1, 1], [-1, -1, -1]]
        self.assertEqual(generate_frequencies(k3), freqs_k3)


if __name__ == '__main__':
    unittest.main()
