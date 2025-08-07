import matplotlib.pyplot as plt

from irispy.tests.helpers import figure_test


@figure_test
def test_simple_plot():
    plt.plot([0, 1])
