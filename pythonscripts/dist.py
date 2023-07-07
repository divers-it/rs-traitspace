from sklearn.metrics.pairwise import nan_euclidean_distances
import pandas as pd

df = pd.read_table('./outputs/diaz_pcf.txt', sep="\t")  

print(df)

df_dist = nan_euclidean_distances(df) # distance between rows of X

df_dist.to_csv('df_dist.csv')