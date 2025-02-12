import numpy as np
import sys
import RNA
import argparse
import subprocess
import pandas as pd
from typing import List
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import seaborn as sns

def calculate_bpp_in_memory(seq):
    """Calculate base pair probabilities (BPP) for a given RNA sequence."""
    print("-> Calculating base-pairing probability...")

    fc = RNA.fold_compound(seq)
    fc.pf()
    basepair_probs = fc.bpp()
    bpp_sums = np.zeros(len(seq))

    for i in range(1, len(basepair_probs)):
        for j in range(i + 1, len(basepair_probs[i])):
            prob = basepair_probs[i][j]
            if prob > 0:
                bpp_sums[i-1] += prob
                bpp_sums[j-1] += prob

    return basepair_probs, bpp_sums

def find_bpp_regions_from_array(bpp_sums, threshold):
    """Identify regions where BPP exceeds the threshold for at least 10 consecutive bases."""
    regions = []
    count = 0
    start = -1
    for i, prob_sum in enumerate(bpp_sums):
        if prob_sum > threshold:
            if count == 0:
                start = i
            count += 1
        else:
            if count >= 10:
                center = start + count // 2
                regions.append(center)
            count = 0

    if count >= 10:
        center = start + count // 2
        regions.append(center)

    return regions

def create_penalty_regions(regions):
    """Generate a list of penalty regions based on identified BPP regions."""
    penalty_regions = []
    for N in regions:
        penalty_regions.extend([
            f"{N}~{N}",
            f"{N}~{N+1}",
            f"{N-1}~{N}",
            f"{N-1}~{N+1}",
            f"{N-2}~{N+2}",
        ])
    return penalty_regions

def get_pairs(s):
    """Extract base pairs from a secondary structure string."""
    pairs = []
    stack = []
    for i, c in enumerate(s):
        if c == '(':
            stack.append(i)
        elif c == ')':
            pairs.append((stack.pop(), i))
    return pairs

def calc_jaccard_folding(fold1, fold2):
    """Compute the Jaccard Index between two RNA secondary structures."""
    assert len(fold1) == len(fold2)

    pairs1 = set(get_pairs(fold1))
    pairs2 = set(get_pairs(fold2))

    return len(pairs1 & pairs2) / len(pairs1 | pairs2)

def calculate_jaccard_for_all_sequences(seq_control: str, sequences: List[str], df: pd.DataFrame) -> List[float]:
    """Calculate the Jaccard Index for a control sequence against all other sequences."""
    jaccard_indices = []
    fold_control = df[df['Sequence'] == seq_control]['Structure'].values[0]
    
    for seq in sequences:
        fold_seq = df[df['Sequence'] == seq]['Structure'].values[0]
        jaccard_index = calc_jaccard_folding(fold_control, fold_seq)
        jaccard_indices.append(round(jaccard_index, 3))
        
    return jaccard_indices

def parse_sample_file_to_df(filename):
    """Parse a LinearDesign sample file and return it as a DataFrame."""
    data = {
        'Name': [],
        'Penalty': [],
        'MFE': [],
        'CAI': [],
        'JI': [],
        'Sequence': [],
        'Structure': [],
    }

    with open(filename, 'r') as file:
        blocks = file.read().strip().split('>')[1:]

        for block in blocks:
            lines = block.strip().split('\n')
            name = lines[0].strip()
            penalty = lines[1].split(': ')[1].strip()
            sequence = lines[2].split(': ')[1].strip()
            structure, mfe = RNA.fold(sequence)
            mfe = round(mfe, 3)
            cai = float(lines[5].split(': ')[1].strip())

            data['Name'].append(name)
            data['Penalty'].append(penalty)
            data['Sequence'].append(sequence)
            data['Structure'].append(structure)
            data['MFE'].append(mfe)
            data['CAI'].append(cai)
            data['JI'].append(-1)

    df = pd.DataFrame(data)
    return df

def parse_and_analyze_results(title, jaccard_threshold):
    """Parse the LinearDesign output file and compute the Jaccard Index for all sequences."""
    print("-> Analyzing samples...")

    df = parse_sample_file_to_df(f"{title}.txt")
    df = df.drop_duplicates(subset='Sequence', keep='first')

    try:
        seq_control = df[df['Name'].str.endswith('with_utr')]['Sequence'].values[0]
    except IndexError:
        print("No control sequence found with 'with_utr' in the 'Name'.")
        seq_control = None

    if seq_control:
        sequences = df['Sequence'].tolist()
        df['JI'] = calculate_jaccard_for_all_sequences(seq_control, sequences, df)
    else:
        df['JI'] = [None] * len(df)
    
    df = df.sort_values(by=['JI'], ascending=[False])
    df = df[df['JI'] <= jaccard_threshold]
    
    df.to_csv(f"{args.output}.tsv", sep='\t', index=False)

def draw_bpp_dotplot(seq_length, output_prefix):
    print("-> Drawing a dotplot...")

    matrix = [[0.0 for _ in range(seq_length)] for _ in range(seq_length)]

    bpp_file = f"{output_prefix}_bpp.txt"
    with open(bpp_file, 'r') as f:
        for line in f.readlines():
            line = line.strip()
            if line == "": 
                continue
            i, j, prob = line.split()
            i, j = int(i) - 1, int(j) - 1
            if 0 <= i < seq_length and 0 <= j < seq_length:
                matrix[j][i] = float(prob)

    matrix_ticks = pd.DataFrame(
        data=matrix,
        columns=range(1, seq_length + 1),
        index=range(1, seq_length + 1)
    )
    mask = np.triu(np.ones_like(matrix, dtype=bool))

    plt.figure(figsize=(50, 50))
    sns.set_theme(style="white")
    heatmap = sns.heatmap(
        matrix_ticks,
        mask=mask,
        cmap="rocket_r",
        vmax=1.0,
        vmin=0.0,
        square=True,
        linewidths=0.5,
        cbar_kws={"shrink": 0.8},
        xticklabels=25,
        yticklabels=25,
    )

    ax = plt.gca()
    ax.tick_params(axis='x', labelsize=20)
    ax.tick_params(axis='y', labelsize=20)

    # Add a black border around the entire heatmap and color bar area
    for _, spine in ax.spines.items():
        spine.set_edgecolor('black')
        spine.set_linewidth(2)

    cbar = heatmap.collections[0].colorbar
    cbar.set_label('Base Pair Probability', fontsize=30, labelpad=20)
    cbar.ax.tick_params(labelsize=30)

    plt.title(f"{output_prefix.split('/')[-1].strip()}: BPP dotplot of original sequence", fontsize=50)
    plt.savefig(f"{output_prefix}_dotplot.png", bbox_inches='tight', pad_inches=2)
    plt.close()


def draw_scatter_plot(output_prefix):
    """Create and save a scatter plot based on CAI, MFE, and Jaccard Index."""
    print("-> Drawing a scatter plot...")
    df = pd.read_csv(f"{output_prefix}.tsv", sep='\t')
    
    plt.figure(figsize=(10, 6))

    colors = ["#000000", "#1E4174", "#6A7BA2", "#FCF6F5"]
    custom_cmap = LinearSegmentedColormap.from_list("custom_cmap", colors)

    scatter = plt.scatter(df['CAI'], df['MFE'], c=df['JI'], cmap=custom_cmap, alpha=0.7, edgecolor='k')
    scatter.set_clim(0.3, 1)

    cbar = plt.colorbar(scatter)
    cbar.set_label('Jaccard Index (JI)')
    
    plt.xlabel('CAI (CDS)')
    plt.ylabel('MFE (kcal/mol)')
    plt.title(f'Samples of {output_prefix}')

    plt.xlim(df['CAI'].min() - 0.01, df['CAI'].max() + 0.01)
    plt.ylim(df['MFE'].min() - 3, df['MFE'].max() + 3)

    plt.grid(True)
    plt.title(f"{output_prefix.split("/")[-1].strip()}", fontsize=15)
    plt.savefig(f"{output_prefix}.png")

def run_linear_design(output_prefix, sequence, utr3, lam, penalty_regions):
    """Run the LinearDesign algorithm with specified penalty regions."""
    print("-> Generating samples...")
    name = output_prefix.split("/")[-1].strip()

    print(f"   ├─ original string")
    subprocess.run(
        f'echo "> {name}_with_utr\nPenalty: -1" > {output_prefix}.txt',
        shell=True
    )
    subprocess.run(
        f'echo "{sequence}" | ./lineardesign -l {lam} -u3 {utr3} >> {output_prefix}.txt',
        shell=True
    )

    for penalty in penalty_regions:
        print(f"   ├─ penalty on {penalty}")
        subprocess.run(
            f'echo "> {name}_{penalty}" >> {output_prefix}.txt',
            shell=True
        )
        subprocess.run(
            f'echo "{sequence}" | ./lineardesign -l {lam} -u3 {utr3} -p {penalty} >> {output_prefix}.txt',
            shell=True
        )
    print("   └─ done.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate RNA base pair probability and analyze LinearDesign results.")
    parser.add_argument("-o", "--output", type=str, default="samples", help="Output file prefix")
    parser.add_argument("-e", "--eval", action='store_true', help="Evaluate original result (BPP, dotplot)")
    parser.add_argument("-t", "--threshold", type=float, default=0.999, help="Threshold for base pair probabilities")
    parser.add_argument("-m", "--manual", type=str, default="", help="Manually provide penalty points")
    parser.add_argument("-l", "--lambda_val", type=float, default=0, help="Lambda value for LinearDesign")
    parser.add_argument("-u3", "--utr3", type=str, default="", help="3' UTR length")
    parser.add_argument("-j", "--jaccard-index", type=float, default=1.0, help="Jaccard Index threshold")
    parser.add_argument("-s", "--scatter", action='store_true', help="Draw scatter plot")
    
    args = parser.parse_args()
    seq = sys.stdin.read().strip()
    
    if not seq:
        print("Error: No RNA sequence provided.")
        sys.exit(1)

    lambda_val = str(args.lambda_val)
    eval_mode = args.eval
    utr3 = str(args.utr3)
    scatter = 1 if args.scatter else 0
    
    if args.manual == "":
        print("-> Running LinearDesign...")
        result = subprocess.run(
            f'echo "{seq}" | ./lineardesign -l {lambda_val} -u3 {utr3}',
            shell=True, capture_output=True, text=True
        )
        sequence_line = result.stdout.split('\n')[-6]
        new_sequence = sequence_line.split(':')[1].strip()

        basepair_probs, bpp_sums = calculate_bpp_in_memory(new_sequence)

        if eval_mode:
            with open(f"{args.output}_bpp.txt", "w") as f:
                for i in range(1, len(basepair_probs)):
                    for j in range(i + 1, len(basepair_probs[i])):
                        prob = basepair_probs[i][j]
                        if prob > 0:
                            f.write(f"{i}\t{j}\t{prob}\n")
            draw_bpp_dotplot(len(new_sequence), args.output)
            sys.exit(0)

        regions = find_bpp_regions_from_array(bpp_sums, args.threshold)
        print(f"-> Penalty points based on BPP (threshold: {args.threshold}): {regions}")
    else:
        regions = list(map(int, args.manual.split(',')))
        print(f"-> Penalty points (manually seleted): {regions}")

    penalty_regions = create_penalty_regions(regions)
    run_linear_design(args.output, seq, utr3, lambda_val, penalty_regions)
    parse_and_analyze_results(args.output, args.jaccard_index)

    if scatter:
        draw_scatter_plot(args.output)
