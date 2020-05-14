import seaborn as sns; sns.set(color_codes=True)
from matplotlib import pyplot as plt

# print(iris["sepal_length"])
# for x in iris: print(x)

species = iris.pop("species")

lut = dict(zip(species.unique(), "rbg"))
row_colors = species.map(lut)
g = sns.clustermap(iris, row_colors=row_colors)
plt.savefig('practice.png',
                dpi=150, figsize=(10, 10))