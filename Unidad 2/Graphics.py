import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from wordcloud import WordCloud

csv_file = "kmer_frequencies.csv"
df = pd.read_csv(csv_file)

sns.set_theme(style="whitegrid")

sorted_df = df.sort_values(by="Frequency", ascending=False)

single_color = "royalblue" 

plt.figure(figsize=(10, 6))
top10 = sorted_df.head(10)
sns.barplot(x="k-mer", y="Frequency", data=top10, color=single_color)
plt.title("Top 10 k-mers mas frecuentes", fontsize=14, fontweight='bold')
plt.xlabel("k-mer", fontsize=12)
plt.ylabel("Frecuencia", fontsize=12)
plt.xticks(rotation=45, fontsize=11)
plt.tight_layout()
plt.savefig("top10_kmers.png", dpi=300)
plt.show()

plt.figure(figsize=(14, 6))
top20 = sorted_df.head(20)
sns.barplot(x="k-mer", y="Frequency", data=top20, color=single_color)
plt.title("Top 20 k-mers mas frecuentes", fontsize=14, fontweight='bold')
plt.xlabel("k-mer", fontsize=12)
plt.ylabel("Frecuencia", fontsize=12)
plt.xticks(rotation=90, fontsize=11)
plt.tight_layout()
plt.savefig("top20_kmers.png", dpi=300)
plt.show()

plt.figure(figsize=(12, 6))
sns.histplot(data=df, x="Frequency", bins=60, color=single_color, edgecolor="black")
plt.yscale("log") 
plt.title("Distribucion de Frecuencias de k-mers (Eje Y Logaritmico)", fontsize=14, fontweight='bold')
plt.xlabel("Frecuencia", fontsize=12)
plt.ylabel("Conteo de k-mers (Escala Logarítmica)", fontsize=12)
plt.tight_layout()
plt.savefig("kmer_histogram.png", dpi=300)
plt.show()

word_freq = dict(zip(df["k-mer"], df["Frequency"]))

wordcloud = WordCloud(
    width=1200,
    height=600,
    background_color='white',
    colormap='plasma' 
).generate_from_frequencies(word_freq)

plt.figure(figsize=(14, 7))
plt.imshow(wordcloud, interpolation='bilinear')
plt.axis("off")
plt.title("Nube de Palabras Genética", fontsize=16, fontweight='bold')
plt.tight_layout()
plt.savefig("kmer_wordcloud.png", dpi=300)
plt.show()

print("\n===== ESTADISTICAS BASICAS =====")
print(f"Total de k-mers unicos: {len(df)}")
print(f"Frecuencia maxima: {df['Frequency'].max()}")
print(f"Frecuencia minima: {df['Frequency'].min()}")
print(f"Frecuencia promedio: {df['Frequency'].mean():.2f}")
print(f"Frecuencia mediana: {df['Frequency'].median():.2f}")

print("\nArchivos generados:")
print("- top10_kmers.png")
print("- top20_kmers.png")
print("- kmer_histogram.png")
print("- kmer_wordcloud.png")