import unittest


class TestVersion(unittest.TestCase):
    def test_version(self):
        """
        Example test
        """
        from oepandas_mae import __version__
        self.assertIsNotNone(__version__)
