import os
import pandas as pd
import matplotlib.pyplot as plt
from adjustText import adjust_text



# Function to generate scatter graph
def create_scatter_graph(title, xlabel, ylabel, y_values, file_name, kinase_avg_AFConf):
    plt.figure(figsize=(10, 6))
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)

    scatter_objects = []
    texts = []

    for kinase_name, y_value in y_values.items():
        avg_AFConf = kinase_avg_AFConf[kinase_name]
        scatter = plt.scatter(avg_AFConf, y_value, marker='x')
        scatter_objects.append(scatter)

        text = plt.text(avg_AFConf, y_value, kinase_name, fontsize=9, ha='right', va='bottom', alpha=0.7)
        texts.append(text)

    x = [kinase_avg_AFConf[kinase_name] for kinase_name in y_values.keys()]
    y = list(y_values.values())
    try:
        adjust_text(texts, x=x, y=y, arrowprops=dict(arrowstyle="-", color='black', lw=0.5), autoalign='xy', only_move={'points':'y'})
    except:
        adjust_text(texts, x=x, y=y, arrowprops=dict(arrowstyle="-", color='black', lw=0.5), autoalign='xy', only_move={'points':'y'})

    plt.tight_layout()
    savefile = os.path.join(folder_path, file_name)
    plt.savefig(savefile)
    plt.show()


# Path to the directory containing kinase files
folder_path = r"/mnt/c/Users/crisd/Desktop/ProteinDesign/results/human_10_gene_min"

# Initialize dictionaries to store kinase data
kinase_avg_alpha = {}
kinase_avg_beta = {}
kinase_avg_AFConf = {}
kinase_avg_score = {}
kinase_avg_nextaa = {}
kinase_mode_score = {}
kinase_median_score = {}
kinase_mode_nextaa = {}
kinase_median_nextaa = {}
kinase_frequence_s={}
kinase_frequence_t={}
kinase_frequence_y={}
kinase_frequence_h={}
kinase_frequence_r={}
kinase_frequence_k={}
data_groups = {
    "in_vivo": {"score": [], "AFConf": [], "alpha": [], "beta": [], "hydturn": [],"NextAA":[],"site_aa":[]},
    "in_vitro": {"score": [], "AFConf": [], "alpha": [], "beta": [], "hydturn": [],"NextAA":[],"site_aa":[]},
    "both": {"score": [], "AFConf": [], "alpha": [], "beta": [], "hydturn": [],"NextAA":[],"site_aa":[]},
    "either": {"score": [], "AFConf": [], "alpha": [], "beta": [], "hydturn": [],"NextAA":[],"site_aa":[]}
}
site_groups = {
    "Y": {"score": [], "AFConf": [], "alpha": [], "beta": [], "NextAA": []},
    "S": {"score": [], "AFConf": [], "alpha": [], "beta": [], "NextAA": []},
    "T": {"score": [], "AFConf": [], "alpha": [], "beta": [], "NextAA": []},
    "R": {"score": [], "AFConf": [], "alpha": [], "beta": [], "NextAA": []},
    "K": {"score": [], "AFConf": [], "alpha": [], "beta": [], "NextAA": []},
    "H": {"score": [], "AFConf": [], "alpha": [], "beta": [], "NextAA": []}
}

# Process each kinase file in the directory
for file_name in os.listdir(folder_path):
    if file_name.endswith("_human_results.txt"):
        kinase_name = file_name.split("_")[0]
        file_path = os.path.join(folder_path, file_name)
        column_names = ["Gene", "Res", "Score",  "alpha", "beta", "iso_B", "alpha3", "alphaI", "hydturn", "AFConf", "NextAA", "orga", "site_aa","sub_acc_id", "in_vivo", "in_vitro"]
        data = pd.read_csv(file_path, sep="\t", names=column_names, skiprows=1, index_col=False)
        data_length =data['Gene'].size
        if data_length  < 50:
            continue

        kinase_avg_alpha[kinase_name] = data["alpha"].mean()
        kinase_avg_beta[kinase_name] = data["beta"].mean()
        
        total_letters = len(data)
        kinase_frequence_s[kinase_name] = data["site_aa"].value_counts().get('S', 0) * 100/data_length 
        kinase_frequence_t[kinase_name] = data["site_aa"].value_counts().get('T', 0) * 100/data_length 
        kinase_frequence_y[kinase_name] = data["site_aa"].value_counts().get('Y', 0) * 100/data_length 
        kinase_frequence_r[kinase_name] = data["site_aa"].value_counts(normalize=True).get('R', 0) * 100
        kinase_frequence_h[kinase_name] = data["site_aa"].value_counts(normalize=True).get('H', 0) * 100
        kinase_frequence_k[kinase_name] = data["site_aa"].value_counts(normalize=True).get('K', 0) * 100

        kinase_mode_score[kinase_name] = data["Score"].mode()[0] if not data["Score"].mode().empty else None
        kinase_median_score[kinase_name] = data["Score"].median()
        kinase_mode_nextaa[kinase_name] = data["NextAA"].mode()[0] if not data["Score"].mode().empty else None
        kinase_median_nextaa[kinase_name] = data["NextAA"].median()
    
        i = 0
        AFsum = 0
        nextaasum=0
        for AFconf in data["AFConf"]:
            if AFconf is not None and 0 < AFconf < 100:
                AFsum += AFconf
                i += 1
        kinase_avg_AFConf[kinase_name] = AFsum / i
        i=0
        for nextaa in data["NextAA"]:
            if nextaa <1500:
                nextaasum += nextaa
                i += 1
        kinase_avg_nextaa[kinase_name] = nextaasum / i
        kinase_avg_score[kinase_name] = data["Score"].mean()
        for index, row in data.iterrows():
            in_vivo = row["in_vivo"]
            in_vitro = row["in_vitro"]
            score = row["Score"]
            AFConf = row["AFConf"]
            alpha = row["alpha"]
            beta = row["beta"]
            hydturn = row["hydturn"]
            NextAA=row["NextAA"]
            site_aa=row["site_aa"]

            if in_vivo == "X" and in_vitro == "X":
                data_groups["both"]["score"].append(score)
                data_groups["both"]["AFConf"].append(kinase_avg_AFConf[kinase_name])
                data_groups["both"]["alpha"].append(alpha)
                data_groups["both"]["beta"].append(beta)
                data_groups["both"]["hydturn"].append(hydturn)  
                data_groups["both"]["NextAA"].append(NextAA)
                data_groups["both"]["site_aa"].append(site_aa)
            elif in_vivo == "X":
                data_groups["in_vivo"]["score"].append(score)
                data_groups["in_vivo"]["AFConf"].append(kinase_avg_AFConf[kinase_name])
                data_groups["in_vivo"]["alpha"].append(alpha)
                data_groups["in_vivo"]["beta"].append(beta)
                data_groups["in_vivo"]["hydturn"].append(hydturn)
                data_groups["in_vivo"]["NextAA"].append(NextAA)
                data_groups["in_vivo"]["site_aa"].append(site_aa)
            elif in_vitro == "X":
                data_groups["in_vitro"]["score"].append(score)
                data_groups["in_vitro"]["AFConf"].append(kinase_avg_AFConf[kinase_name])
                data_groups["in_vitro"]["alpha"].append(alpha)
                data_groups["in_vitro"]["beta"].append(beta)
                data_groups["in_vitro"]["hydturn"].append(hydturn)
                data_groups["in_vitro"]["NextAA"].append(NextAA)
                data_groups["in_vitro"]["site_aa"].append(site_aa)
            if in_vivo == "X" or in_vitro == "X":
                data_groups["either"]["score"].append(score)
                data_groups["either"]["AFConf"].append(kinase_avg_AFConf[kinase_name])
                data_groups["either"]["alpha"].append(alpha)
                data_groups["either"]["beta"].append(beta)
                data_groups["either"]["hydturn"].append(hydturn)
                data_groups["either"]["NextAA"].append(NextAA)
                data_groups["either"]["site_aa"].append(site_aa)
                site_groups[site_aa]["score"].append(score)
                site_groups[site_aa]["AFConf"].append(AFConf)
                site_groups[site_aa]["alpha"].append(alpha)
                site_groups[site_aa]["beta"].append(beta)
                site_groups[site_aa]["NextAA"].append(NextAA)

# Create a bar plot comparing average alpha and beta for each kinase
kinase_names = list(kinase_avg_alpha.keys())
average_alphas = list(kinase_avg_alpha.values())
average_betas = list(kinase_avg_beta.values())

x = list(range(len(kinase_names)))
width = 0.4

fig, ax = plt.subplots(figsize=(12, 6))
rects1 = ax.bar(x, average_alphas, width, label='Average Alpha')
rects2 = ax.bar([i + width for i in x], average_betas, width, label='Average Beta')

ax.set_ylabel('Averages')
ax.set_title('Average Alpha and Beta for Kinases')
ax.set_xticks([i + width/2 for i in x])
ax.set_xticklabels(kinase_names, rotation=90)
ax.legend()

plt.tight_layout()
savefile=os.path.join(folder_path,"average_alpha_beta_comparison.png")
plt.savefig(savefile)
plt.show()

# Create a scatter plot comparing average AFConf vs. average Score for each kinase
plt.figure(figsize=(10, 6))
plt.title("Average AFConf vs. Average Score for Kinases")
plt.xlabel("Average AFConf")
plt.ylabel("Average Score")

scatter_objects = []

for kinase_name, avg_score in kinase_avg_score.items():
    avg_AFConf = kinase_avg_AFConf[kinase_name]
    scatter = plt.scatter(avg_AFConf, avg_score, marker='x')
    scatter_objects.append(scatter)

# Create a list to store Text objects for annotations
texts = []
for scatter, kinase_name in zip(scatter_objects, kinase_names):
    avg_AFConf = kinase_avg_AFConf[kinase_name]
    avg_score = kinase_avg_score[kinase_name]
    text = plt.text(avg_AFConf, avg_score, kinase_name, fontsize=9, ha='right', va='bottom', alpha=0.7)
    texts.append(text)

# Use adjust_text to automatically adjust the positions of the labels
x = [kinase_avg_AFConf[kinase_name] for kinase_name in kinase_names]
y = [kinase_avg_score[kinase_name] for kinase_name in kinase_names]

adjust_text(texts, x=x, y=y, arrowprops=dict(arrowstyle="-", color='black', lw=0.5), autoalign='xy', only_move={'points':'y', 'text':'y'})

plt.tight_layout()

savefile = os.path.join(folder_path, "average_AFConf_vs_Score.png")
plt.savefig(savefile)
plt.show()


# Create a bar plot comparing different data groups
group_names = ["in_vivo", "in_vitro", "both","either"]
data_names = ["score", "AFConf", "alpha", "beta", "hydturn"]

fig, axes = plt.subplots(nrows=len(data_names), ncols=1, figsize=(10, 12), sharex=True)

for i, data_name in enumerate(data_names):
    ax = axes[i]
    data_values = [data_groups[group][data_name] for group in group_names]
    ax.bar(group_names, [sum(data) / len(data) for data in data_values])
    ax.set_ylabel(f'Average {data_name.capitalize()}')
    ax.set_title(f'Average {data_name.capitalize()} for Different Groups')

plt.xlabel('Data Group')
plt.xticks(ticks=list(range(len(group_names))), labels=group_names)  # Set group names as x-axis ticks
plt.tight_layout()
savefile = os.path.join(folder_path, "average_data_groups_comparison.png")
plt.savefig(savefile)
plt.show()

# Create scatter plots for each measure
create_scatter_graph("Average AFConf vs. Average nextAA for Kinases", "Average AFConf", "Average nextAA", kinase_avg_nextaa, "average_AFConf_vs_AverageNextAA.png", kinase_avg_AFConf)
create_scatter_graph("Average AFConf vs. Mode Score for Kinases", "Average AFConf", "Mode Score", kinase_mode_score, "average_AFConf_vs_ModeScore.png", kinase_avg_AFConf)
create_scatter_graph("Average AFConf vs. Median Score for Kinases", "Average AFConf", "Median Score", kinase_median_score, "average_AFConf_vs_MedianScore.png", kinase_avg_AFConf)
create_scatter_graph("Average AFConf vs. Mode nextAA for Kinases", "Average AFConf", "Mode nextAA", kinase_mode_nextaa, "average_AFConf_vs_ModeNextAA.png", kinase_avg_AFConf)
create_scatter_graph("Average AFConf vs. Median nextAA for Kinases", "Average AFConf", "Median nextAA", kinase_median_nextaa, "average_AFConf_vs_MedianNextAA.png", kinase_avg_AFConf)

kinase_names = list(kinase_avg_alpha.keys())
frequencies_s = list(kinase_frequence_s.values())
frequencies_t = list(kinase_frequence_t.values())
frequencies_k = list(kinase_frequence_k.values())
frequencies_r = list(kinase_frequence_r.values())
frequencies_y = list(kinase_frequence_y.values())
frequencies_h = list(kinase_frequence_h.values())

x = list(range(len(kinase_names)))
width = 0.23

fig, ax = plt.subplots(figsize=(18, 8))

rects_s = ax.bar([i for i in x], frequencies_s, width, label='S')
rects_t = ax.bar([i - width for i in x], frequencies_t, width, label='T')
rects_y = ax.bar([i + width for i in x], frequencies_y, width, label='Y')


ax.set_ylabel('Letter Frequencies (%)')
ax.set_title('Letter Frequencies for Kinases')
ax.set_xticks(x)
ax.set_xticklabels(kinase_names, rotation=90)
ax.legend()

plt.tight_layout()
savefile = os.path.join(folder_path, "letter_frequencies_comparison.png")
plt.savefig(savefile)
plt.show()

# Create a bar plot comparing different measures for different site groups
site_group_names = ["Y", "S", "T"]
data_names = ["AFConf", "score", "alpha", "beta", "NextAA"]

fig, axes = plt.subplots(nrows=len(data_names), ncols=1, figsize=(10, 15), sharex=True)

for i, data_name in enumerate(data_names):
    ax = axes[i]
    data_values = []
    for site_group in site_group_names:
        data_values.append(site_groups[site_group][data_name])
    
    ax.bar(site_group_names, [sum(data) / len(data) for data in data_values])
    ax.set_ylabel(f'Average {data_name.capitalize()}')
    ax.set_title(f'Average {data_name.capitalize()} for Different Site Groups')

plt.tight_layout()
savefile = os.path.join(folder_path, "average_site_groups_comparison.png")
plt.savefig(savefile)
plt.show()