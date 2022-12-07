import pandas as pd
import matplotlib.pyplot as plt

df=pd.read_csv('Lax.csv')
df.plot()
plt.title('Lax')
plt.savefig('Lax.png')

df=pd.read_csv('L_W.csv')
df.plot()
plt.title('L-W')
plt.savefig('L_W.png')

df=pd.read_csv('Mac.csv')
df.plot()
plt.title('MacC')
plt.savefig('Mac.png')

df=pd.read_csv('FVS.csv')
df.plot()
plt.title('FVS')
plt.savefig('FVS.png')

df=pd.read_csv('Roe.csv')
df.plot()
plt.title('Roe')
plt.savefig('Roe.png')

#plt.show()