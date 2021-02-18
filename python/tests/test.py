import os
import pathlib
import unittest

import numpy as np

from libsonatareport import (SpikeReader, SpikePopulation,
                       SomaReportReader, SomaReportPopulation,
                       ElementReportReader, ElementReportPopulation,
                       )


PATH = os.path.dirname(os.path.realpath(__file__))
PATH = os.path.join(PATH, '../../tests/data')


class TestSpikePopulation(unittest.TestCase):
    def setUp(self):
        path = os.path.join(PATH, "spikes.h5")
        self.test_obj = SpikeReader(path)

    def test_get_all_populations(self):
        self.assertEqual(self.test_obj.get_population_names(), ['All', 'empty', 'spikes1', 'spikes2'])

    def test_get_population(self):
        self.assertTrue(isinstance(self.test_obj['spikes1'], SpikePopulation))

    def test_get_inexistant_population(self):
        self.assertRaises(RuntimeError, self.test_obj.__getitem__, 'foobar')

    def test_get_spikes_from_population(self):
        self.assertEqual(self.test_obj['All'].get(), [(5, 0.1), (2, 0.2), (3, 0.3), (2, 0.7), (3, 1.3)])
        self.assertEqual(self.test_obj['All'].get(tstart=0.2, tstop=1.0), [(2, 0.2), (3, 0.3), (2, 0.7)])
        self.assertEqual(self.test_obj['spikes2'].get(tstart=0.2, tstop=1.0), [(3, 0.3), (2, 0.2), (2, 0.7)])
        self.assertEqual(self.test_obj['spikes1'].get((3,)), [(3, 0.3), (3, 1.3)])
        self.assertEqual(self.test_obj['spikes2'].get((3,)), [(3, 0.3), (3, 1.3)])
        self.assertEqual(self.test_obj['spikes2'].get((10,)), [])
        self.assertEqual(self.test_obj['spikes2'].get((2,), 0., 0.5), [(2, 0.2)])
        self.assertEqual(self.test_obj['spikes1'].get((2, 5)), [(2, 0.2), (2, 0.7), (5, 0.1)])
        self.assertEqual(self.test_obj['spikes2'].get((2, 5)), [(5, 0.1), (2, 0.2), (2, 0.7)])
        self.assertEqual(self.test_obj['All'].sorting, "by_time")
        self.assertEqual(self.test_obj['spikes1'].sorting, "by_id")
        self.assertEqual(self.test_obj['spikes2'].sorting, "none")
        self.assertEqual(self.test_obj['empty'].get(), [])

        self.assertEqual(len(self.test_obj['All'].get(node_ids=[])), 0)

class TestSomaReportPopulation(unittest.TestCase):
    def setUp(self):
        path = os.path.join(PATH, "somas.h5")
        self.test_obj = SomaReportReader(path)

    def test_get_all_population(self):
        self.assertEqual(self.test_obj.get_population_names(), ['All', 'soma1', 'soma2'])

    def test_get_population(self):
        self.assertTrue(isinstance(self.test_obj['All'], SomaReportPopulation))

    def test_get_inexistant_population(self):
        self.assertRaises(RuntimeError, self.test_obj.__getitem__, 'foobar')

    def test_get_reports_from_population(self):
        self.assertEqual(self.test_obj['All'].times, (0., 1., 0.1))
        self.assertEqual(self.test_obj['All'].time_units, 'ms')
        self.assertEqual(self.test_obj['All'].data_units, 'mV')
        self.assertFalse(self.test_obj['All'].sorted)
        self.assertEqual(len(self.test_obj['All'].get().ids), 20)  # Number of nodes
        self.assertEqual(len(self.test_obj['All'].get().times), 10)  # number of times
        self.assertEqual(len(self.test_obj['All'].get().data), 10)  # should be the same

        sel = self.test_obj['All'].get(node_ids=[13, 14], tstart=0.8, tstop=1.0)
        self.assertEqual(len(sel.times), 2)  # Number of timestamp (0.8 and 0.9)
        self.assertEqual(list(sel.ids), [13, 14])
        np.testing.assert_allclose(sel.data, [[13.8, 14.8], [13.9, 14.9]])

        sel_all = self.test_obj['All'].get()
        self.assertEqual(sel_all.ids, [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20])

        sel_empty = self.test_obj['All'].get(node_ids=[])
        np.testing.assert_allclose(sel_empty.data, np.empty(shape=(0, 0)))


class TestElementReportPopulation(unittest.TestCase):
    def setUp(self):
        path = os.path.join(PATH, "elements.h5")
        self.test_obj = ElementReportReader(path)

    def test_get_all_population(self):
        self.assertEqual(self.test_obj.get_population_names(), ['All', 'element1', 'element42'])

    def test_get_population(self):
        self.assertTrue(isinstance(self.test_obj['All'], ElementReportPopulation))

    def test_get_inexistant_population(self):
        self.assertRaises(RuntimeError, self.test_obj.__getitem__, 'foobar')

    def test_get_node_ids(self):
        self.assertEqual(self.test_obj['All'].get_node_ids(), [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20])

    def test_get_reports_from_population(self):
        self.assertEqual(self.test_obj['All'].times, (0., 4., 0.2))
        # check following calls succeed (no memory destroyed)
        self.assertEqual(self.test_obj['All'].times, (0., 4., 0.2))

        self.assertEqual(self.test_obj['All'].time_units, 'ms')
        self.assertEqual(self.test_obj['All'].data_units, 'mV')
        self.assertTrue(self.test_obj['All'].sorted)
        self.assertEqual(len(self.test_obj['All'].get(tstride=2).data), 10)  # Number of times in this range
        self.assertEqual(len(self.test_obj['All'].get(tstride=2).times), 10)  # Should be the same
        self.assertEqual(len(self.test_obj['All'].get().ids), 100)
        sel = self.test_obj['All'].get(node_ids=[13, 14], tstart=0.8, tstop=1.2)
        keys = list(sel.ids)
        keys.sort()
        self.assertEqual(keys, [(13, 30), (13, 30), (13, 31), (13, 31), (13, 32), (14, 32), (14, 33), (14, 33), (14, 34), (14, 34)])

        self.assertEqual(len(self.test_obj['All'].get(node_ids=[]).data), 0)
        self.assertEqual(len(self.test_obj['All'].get(node_ids=[]).times), 0)
        self.assertEqual(len(self.test_obj['All'].get(node_ids=[]).ids), 0)

        self.assertEqual(len(sel.times), 3)  # Number of timestamp (0.8, 1.0 and 1.2)
        with self.assertRaises(SonataError):
            self.test_obj['All'].get(tstart=5.)  # tstart out of range

        # tstart should be <= tstop
        np.testing.assert_allclose(self.test_obj['All'].get(node_ids=[1, 2], tstart=3., tstop=3.).data[0],
                                   [150.0, 150.1, 150.2, 150.3, 150.4, 150.5, 150.6, 150.7, 150.8, 150.9])
        # check following calls succeed (no memory destroyed)
        np.testing.assert_allclose(self.test_obj['All'].get(node_ids=[1, 2], tstart=3., tstop=3.).data[0],
                                   [150.0, 150.1, 150.2, 150.3, 150.4, 150.5, 150.6, 150.7, 150.8, 150.9])
        np.testing.assert_allclose(self.test_obj['All'].get(node_ids=[3, 4], tstart=0.2, tstop=0.4).data[0],
                                   [11.0, 11.1, 11.2, 11.3, 11.4, 11.5, 11.6, 11.7, 11.8, 11.9], 1e-6, 0)
        np.testing.assert_allclose(self.test_obj['All'].get(node_ids=[3, 4], tstride=4).data[2],
                                   [81.0, 81.1, 81.2, 81.3, 81.4, 81.5, 81.6, 81.7, 81.8, 81.9], 1e-6, 0)


def test_path_ctor():
    #  make sure constructors that take file paths can use pathlib.Path
    path = pathlib.Path(PATH)

    SpikeReader(path / 'spikes.h5')
    SomaReportReader(path / 'somas.h5')
    ElementReportReader(path / 'elements.h5')


if __name__ == '__main__':
    unittest.main()
