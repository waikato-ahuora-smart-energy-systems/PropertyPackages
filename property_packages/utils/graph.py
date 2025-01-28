import sys
import os
import json
from block_graph import plot_graph, plot_change

def main():
    # Define the path to the /graphs folder
    graphs_folder = './graphs'  # Adjust if necessary

    # if len(sys.argv) > 1:
    #     if sys.argv[0] == "avg":
        
    #     elif sys.argv[0] == "max":

    # Check if the directory exists
    if not os.path.exists(graphs_folder):
        print(f"The directory {graphs_folder} does not exist.")
        sys.exit(1)

    # Find all .json files in the /graphs folder
    json_files = [os.path.join(graphs_folder, f) for f in os.listdir(graphs_folder) if f.endswith('.json')]

    # Check if there are no .json files
    if not json_files:
        print(f"No .json files found in {graphs_folder}.")
        sys.exit(1)

    # Process each .json file
    for json_file in json_files:
        plot_graph(json_file)
    
    # plot_change(json_files)

if __name__ == "__main__":
    main()