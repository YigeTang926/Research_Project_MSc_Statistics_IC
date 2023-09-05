import os
import pandas as pd
import matplotlib.pyplot as plt
plt.style.use('seaborn')

# Get the path
args = {}
args['cwd'] = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir, os.pardir))
args['raw_dir'] = os.path.join(args['cwd'], 'data', 'raw')
args['derived_dir'] = os.path.join(args['cwd'], 'data', 'derived', 'czech', 'gnar_garch_AIC')
args['figures_dir'] = os.path.join(args['cwd'], 'figures', 'czech')

# Read the data
wind_AIC_1 = pd.read_csv(os.path.join(args['derived_dir'], 'wind_AIC_1.csv'),index_col=0)
wind_AIC_2 = pd.read_csv(os.path.join(args['derived_dir'], 'wind_AIC_2.csv'),index_col=0)
wind_AIC_beta_1 = pd.read_csv(os.path.join(args['derived_dir'], 'wind_AIC_beta_1.csv'),index_col=0)
wind_AIC_beta_2 = pd.read_csv(os.path.join(args['derived_dir'], 'wind_AIC_beta_2.csv'),index_col=0)
solar_AIC_1 = pd.read_csv(os.path.join(args['derived_dir'], 'solar_AIC_1.csv'),index_col=0)
solar_AIC_2 = pd.read_csv(os.path.join(args['derived_dir'], 'solar_AIC_2.csv'),index_col=0)
solar_AIC_beta_1 = pd.read_csv(os.path.join(args['derived_dir'], 'solar_AIC_beta_1.csv'),index_col=0)
solar_AIC_beta_2 = pd.read_csv(os.path.join(args['derived_dir'], 'solar_AIC_beta_2.csv'),index_col=0)

# Contour plot
fig, ax = plt.subplots(4,2,figsize=(8,12))
contour1 = ax[0,0].contourf(wind_AIC_1.T,cmap='viridis')
ax[0,0].set_title('$\mathcal{L}_2$ (wind)')
ax[0,0].set_ylabel('$s_1=\cdots=s_p=s$')
ax[0,0].set_yticks(range(7))
ax[0,0].set_xlabel('$p$')
ax[0,0].set_xticks(range(10))
ax[0,0].set_xticklabels(range(1,11))
fig.colorbar(contour1, ax=ax[0,0])
contour2 = ax[0,1].contourf(wind_AIC_2.T,cmap='viridis')
ax[0,1].set_title('$\mathcal{L}_3$ (wind)')
ax[0,1].set_ylabel('$s_1=\cdots=s_p=s$')
ax[0,1].set_yticks(range(7))
ax[0,1].set_xlabel('$p$')
ax[0,1].set_xticks(range(10))
ax[0,1].set_xticklabels(range(1,11))
fig.colorbar(contour2, ax=ax[0,1])
contour3 = ax[1,0].contourf(solar_AIC_1.T,cmap='viridis')
ax[1,0].set_title('$\mathcal{L}_2$ (solar)')
ax[1,0].set_ylabel('$s_1=\cdots=s_p=s$')
ax[1,0].set_yticks(range(7))
ax[1,0].set_xlabel('$p$')
ax[1,0].set_xticks(range(10))
ax[1,0].set_xticklabels(range(1,11))
fig.colorbar(contour3, ax=ax[1,0])
contour4 = ax[1,1].contourf(solar_AIC_2.T,cmap='viridis')
ax[1,1].set_title('$\mathcal{L}_3$ (solar)')
ax[1,1].set_ylabel('$s_1=\cdots=s_p=s$')
ax[1,1].set_yticks(range(7))
ax[1,1].set_xlabel('$p$')
ax[1,1].set_xticks(range(10))
ax[1,1].set_xticklabels(range(1,11))
fig.colorbar(contour4, ax=ax[1,1])
contour5 = ax[2,0].contourf(wind_AIC_beta_1.T,cmap='viridis')
ax[2,0].set_title('$\mathcal{L}_2$ (wind $p=2$)')
ax[2,0].set_ylabel('$s_2$')
ax[2,0].set_yticks(range(7))
ax[2,0].set_xlabel('$s_1$')
ax[2,0].set_xticks(range(7))
fig.colorbar(contour5, ax=ax[2,0])
contour6 = ax[2,1].contourf(wind_AIC_beta_2.T,cmap='viridis')
ax[2,1].set_title('$\mathcal{L}_3$ (wind $p=2$)')
ax[2,1].set_ylabel('$s_2$')
ax[2,1].set_yticks(range(7))
ax[2,1].set_xlabel('$s_1$')
ax[2,1].set_xticks(range(7))
fig.colorbar(contour6, ax=ax[2,1])
contour7 = ax[3,0].contourf(solar_AIC_beta_1.T,cmap='viridis')
ax[3,0].set_title('$\mathcal{L}_2$ (solar $p=2$)')
ax[3,0].set_ylabel('$s_2$')
ax[3,0].set_yticks(range(7))
ax[3,0].set_xlabel('$s_1$')
ax[3,0].set_xticks(range(7))
fig.colorbar(contour7, ax=ax[3,0])
contour8 = ax[3,1].contourf(solar_AIC_beta_2.T,cmap='viridis')
ax[3,1].set_title('$\mathcal{L}_3$ (solar $p=2$)')
ax[3,1].set_ylabel('$s_2$')
ax[3,1].set_yticks(range(7))
ax[3,1].set_xlabel('$s_1$')
ax[3,1].set_xticks(range(7))
fig.colorbar(contour8, ax=ax[3,1])
plt.tight_layout()
plt.savefig(os.path.join(args['figures_dir'], 'gnar_garch_AIC_czech.pdf'))