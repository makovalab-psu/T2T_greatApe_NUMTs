### Saswat's code for upset plots

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
import upsetplot as upsplt

print('# Loaded libaries')

def saswat_generate_upset_plot(
        input_file_path,
        output_file_path,
        kmer_dissim_threshold,
        cardinality=False,
        min_subset_size = 1, # Set the minimum subset size
        max_subset_rank = 1000,
        write=False
    ):
    print(f"[INFO] Reading input file: {input_file_path}")
    # Load and process the presence/absence matrix.
    presAbs_df = pd.read_table(input_file_path).drop('NUMT_ID', axis=1).rename(columns={'Sorang': "S. orangutan", 'Borang': "B. orangutan"})
    print(f"[INFO] Input matrix shape: {presAbs_df.shape}")
    print(f"[INFO] kmer_dissim_threshold={kmer_dissim_threshold}, cardinality={cardinality}, min_subset_size={min_subset_size}, max_subset_rank={max_subset_rank}, write={write}")

    presAbsMatrixNormUpsetDict = {} #empty dictionary to store the indices of the presence/absence matrix.
    for column in presAbs_df.columns: #iterate through the columns.
        indices = [i for i, value in enumerate(presAbs_df[column]) if value == 1] #store the indices of the column where the value is 1.
        presAbsMatrixNormUpsetDict[column] = indices #store the indices in the dictionary.
        print(f"[INFO] Column '{column}': {len(indices)} present entries")

    # Import data into the upset plot.
    upsetData = upsplt.from_contents(presAbsMatrixNormUpsetDict) 
    print(f"[INFO] Raw upset data size: {len(upsetData)}")

    # Filter out subsets that are smaller than the minimum subset size.
    upsetData = upsetData[upsetData >= min_subset_size]
    print(f"[INFO] Filtered upset data size: {len(upsetData)}")

    # Plot with cardinality.
    if cardinality:
        # The real upset plot, however, hidden by making intersection_plot_elements = 0
        upset = upsplt.UpSet(upsetData, sort_by="cardinality", sort_categories_by="-input", facecolor=CBPalette["Black"],
                             max_subset_rank = max_subset_rank, show_counts=False, totals_plot_elements=4, intersection_plot_elements=10)
        print("[INFO] Created UpSet plot with cardinality sorting")
    # Plot by group/clade order.
    else:
        # The real upset plot, however, hidden by making intersection_plot_elements = 0
        upset = upsplt.UpSet(upsetData, sort_categories_by="-input", facecolor=CBPalette["Black"], min_subset_size=min_subset_size,
                            max_subset_rank = max_subset_rank, show_counts=False, totals_plot_elements=4, intersection_plot_elements=10)
        print("[INFO] Created UpSet plot with category ordering")
    
    # Optional color customization (excluding Siamang).
    colorhex = ["#004d71", "#8759a1", "#f75a78", "#ffa600"]
    upset.style_subsets(facecolor=colorhex[0], linewidth=1, present={"Bonobo", "Chimpanzee", "Human", "Gorilla", "S. orangutan", "B. orangutan"})  # BCHGO
    upset.style_subsets(facecolor=colorhex[0], linewidth=1, absent={"S. orangutan", "B. orangutan"}, present={"Bonobo", "Chimpanzee", "Human", "Gorilla"})  # BCHG
    upset.style_subsets(facecolor=colorhex[0], linewidth=1, absent={"S. orangutan", "B. orangutan", "Gorilla"}, present={"Bonobo", "Chimpanzee", "Human"})  # BCH
    upset.style_subsets(facecolor=colorhex[0], linewidth=1, absent={"S. orangutan", "B. orangutan", "Gorilla", "Human"}, present={"Bonobo", "Chimpanzee"})  # BC
    upset.style_subsets(facecolor=colorhex[0], linewidth=1, absent={"Gorilla", "Human", "Bonobo", "Chimpanzee"}, present={"S. orangutan", "B. orangutan"})  # O

    # # Optional color customization (with Siamang).
    # colorhex = ["#004d71", "#8759a1", "#f75a78", "#ffa600"]
    # upset.style_subsets(facecolor=colorhex[0], linewidth=1, absent={"Siamang"}, present={"Bonobo", "Chimpanzee", "Human", "Gorilla", "S. orangutan", "B. orangutan"})  # BCHGO
    # upset.style_subsets(facecolor=colorhex[0], linewidth=1, absent={"Siamang", "S. orangutan", "B. orangutan"}, present={"Bonobo", "Chimpanzee", "Human", "Gorilla"})  # BCHG
    # upset.style_subsets(facecolor=colorhex[0], linewidth=1, absent={"Siamang", "S. orangutan", "B. orangutan", "Gorilla"}, present={"Bonobo", "Chimpanzee", "Human"})  # BCH
    # upset.style_subsets(facecolor=colorhex[0], linewidth=1, absent={"Siamang", "S. orangutan", "B. orangutan", "Gorilla", "Human"}, present={"Bonobo", "Chimpanzee"})  # BC
    # upset.style_subsets(facecolor=colorhex[0], linewidth=1, absent={"Siamang", "Gorilla", "Human", "Bonobo", "Chimpanzee"}, present={"S. orangutan", "B. orangutan"})  # O

    print("[INFO] Rendering plot")
    plot_result = upset.plot() #plot the result.

    # Set the font size of the y-axis tick labels.
    plot_result["matrix"].tick_params(axis='y', labelsize=12)

    # The totals horizontal section.
    plot_result["totals"].xaxis.set_major_formatter(FuncFormatter(comma_formatter))
    plot_result["totals"].set_xlabel("\nTotal # of\nNUMTs in species", fontsize=15, x=0.6) #set the x label
    plot_result["totals"].tick_params(axis='x', labelsize=9) #set the x tick labels

    plt.title(f"Number of NUMTs shared among species\n",fontsize=15)
    plt.ylabel("# of NUMTs\n", fontsize=18) 
    plt.grid(alpha=0.5, linestyle="--")
    print(f"[INFO] Saving plot to: {output_file_path}")
    plt.savefig(output_file_path, format=output_file_path[-3:], transparent=True) #save the plot

    png_output_file_path = output_file_path.rsplit('.', 1)[0] + '.png'
    print(f"[INFO] Saving plot to: {png_output_file_path}")
    plt.savefig(png_output_file_path, format='png', transparent=False)
    

    print("[INFO] Plot complete")
    # plt.show()

print("# Defined main function.")

# Define the formatter function.
def comma_formatter(x, _pos):
    return '{:,.0f}'.format(x)

print("# Defined comma function.")

# Colorblind friendly palette.
CBPalette = {
    "Black": "#000000",
    "Orange": "#F4A637",
    "Light blue": "#B6DBFF",
    "Vermilion": "#DB5829",
    "Mid blue": "#7BB0DF",
    "Maroon": "#894B45",
    "Dark blue": "#1964B0",
    "Light purple": "#D2BBD7",
    "Light teal": "#00C992",
    "Purple": "#AE75A2",
    "Teal": "#008A69",
    "Dark purple": "#882D71",
    "Dark teal": "#386350",
    "Grey": "#DEDEDE",
    "Yellow": "#E9DC6D"
}

print("# Defined color palette.")

#######
# Run:
#######

flanks=500
kmer_dissim_threshold=0.5

for species in [ 'GorGor', 'PanPan', 'PanTro', 'PonAbe', 'PonPyg', 'hprc' ]:
    print(f"[INFO] Starting species: {species}")
    file_in = f"combined_matrix.t2t_{species}.tsv"

    file_out = f"combined_matrix.t2t_{species}.cardinality.pdf"
    saswat_generate_upset_plot(
        input_file_path=file_in,
        output_file_path=file_out,
        kmer_dissim_threshold=kmer_dissim_threshold,
        min_subset_size=1,
        max_subset_rank = 30,
        cardinality=True,
        write=True
    )
    print(f"[INFO] Finished species: {species}")

    # file_out = f"results_{flanks}bp/upset_plot_flanks{flanks}_{kmer_dissim_threshold}.pdf"
    # saswat_generate_upset_plot(
    #     input_file_path=file_in,
    #     output_file_path=file_out,
    #     kmer_dissim_threshold=kmer_dissim_threshold,
    #     min_subset_size=1,
    #     max_subset_rank=1000,
    #     cardinality=False,
    #     write=True
    # )


    # file_out = f"results_{flanks}bp/upset_plot_flanks{flanks}_{kmer_dissim_threshold}_min5.pdf"
    # saswat_generate_upset_plot(
    #     input_file_path=file_in,
    #     output_file_path=file_out,
    #     kmer_dissim_threshold=kmer_dissim_threshold,
    #     min_subset_size=5,
    #     max_subset_rank=1000,
    #     cardinality=False,
    #     write=True
    # )




