# test/test_qplot.py
import unittest
from src.qplot import QPlot

class TestQPlot(unittest.TestCase):

    def setUp(self):
        self.qplot = QPlot()

        # Set Equilibrium-related inputs
        self.qplot.num_eq_entry.set('2')
        for i in range(2):
            self.qplot.eq_pulse_entries[i].set('99776')
            self.qplot.dda_entries[i].set('eftp')
            self.qplot.uid_entries[i].set('jetppf')
            self.qplot.seq_entries[i].set('0')

        # Set TRANSP-related inputs
        self.qplot.num_transp_entry.set('2')
        for i in range(2):
            self.qplot.transp_pulse_entries[i].set('99776')
            self.qplot.runid_entries[i].set('M06')

        # Set other inputs
        self.qplot.time_entry.set('47.4,49.0')
        self.qplot.plot_type_var.set('Profile')

    def test_plot_q_profiles(self):
        # Act
        self.qplot.plot_q_profiles()

        # Assert
        # Here we should assert some conditions that we expect to be true
        # after calling the plot_q_profiles method. This will depend on the
        # specifics of your code.

if __name__ == '__main__':
    unittest.main()
