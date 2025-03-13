import numpy as np
from matplotlib import pyplot as plt
import json

def plot_temperatures():
    # Load the dictionary from the file
    with open("Results/temperatures.json", 'r') as file:
    # with open("probabilities.json", 'r') as file:
        data = json.load(file)

    K = list(data.keys())
    temperatures = list(data.values())

    # Define the sliding window size
    window_size = 3

    # Calculate moving average using the sliding window
    smoothed_temperatures = np.convolve(temperatures, np.ones(window_size) / window_size, mode='valid')

    # Create a discretized X-axis
    discretized_K = np.arange(window_size // 2, len(K) - window_size // 2)

    # Plot the smoothed data
    plt.plot(discretized_K, smoothed_temperatures)
    plt.xlabel('k')
    plt.ylabel('Average Temperature')
    plt.show()

def plot_probabilities():
    # Load the dictionary from the file
    with open("Results/probabilities.json", 'r') as file:
        data = json.load(file)

    K = list(data.keys())
    values = list(data.values())

    # Define the sliding window size
    window_size = 3

    # Calculate moving averages using the sliding window
    smoothed_values = np.convolve(values, np.ones(window_size) / window_size, mode='valid')

    # Create discretized X-axis
    discretized_K = np.arange(window_size // 2, len(K) - window_size // 2)

    # Plot the smoothed data
    plt.plot(discretized_K, smoothed_values)
    plt.xlabel('k')
    plt.ylabel('Average Probability')
    plt.show()

def plot_new_Generated_Mol():
    # Read the JSON file
    with open("Results/dic_news.json", 'r') as file:
        data = json.load(file)

    # Retrieve the keys from the JSON
    keys = list(data.keys())
    max_ite = max(keys)

    for x in range(int(max_ite) + 1):
        if str(x) not in data:
            data[x] = 0

    # Define the window size
    window_size = 5

    # Example data (replace with your actual data)
    iterations = list(data.keys())
    values = list(data.values())

    x_list = []
    y_list = []
    cumulative_sum = 0

    for i in range(len(iterations)):
        cumulative_sum += values[i]
        if i % window_size == 0:
            x_list.append(i)
            y_list.append(cumulative_sum)

    # Plot the smoothed data
    plt.plot(x_list, y_list)
    plt.xlabel('Iteration')
    plt.ylabel('Accumulated quantity of new molecules generated')

    plt.show()

    print(f"Total: {sum(values)}")

def plot_plot_new_Generated_Mol_Criteria():
    with open("Results/dic_news_ok.json", 'r') as file:
        data = json.load(file)

        iterations = list(data.keys())
        accepted_cumulative = list(data.values())

        plt.bar(iterations, accepted_cumulative)
        plt.xlabel('Iteration')
        plt.ylabel('Number of molecules')
        plt.title('Number of New Molecules That Meet the Criteria')
        plt.show()