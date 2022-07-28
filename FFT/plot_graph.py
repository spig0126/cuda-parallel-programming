import matplotlib.pyplot as plt
import pandas as pd

def plot_graph(filepath, format, graph_num, title, x_label, y_label):
    data = pd.read_csv(filepath, sep='\s+', header=None)    #\s = whitespace, \s+ = more than one whitespaces
    data = pd.DataFrame(data)

    x = data[0]
    y = data[1]
    # plt.subplot(graph_num)
    plt.plot(x, y, format)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.title(title)
    plt.show()


if __name__ == '__main__':
    plot_graph("original_signal.txt", "b-", 311, "original_signal", "time", "amplitude")
    plot_graph("idft_signal.txt", "r-", 312, "idft_signal", "time", "amplitude")
    plot_graph("dft_frequencies.txt", "bo", 313, "dft_frequencies", "frequency", "amplitude")


