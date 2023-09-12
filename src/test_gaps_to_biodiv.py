
import torch
import pyro
import pyro.distributions as dist
import geopandas as gpd
import pandas as pd
from pyro.infer import MCMC, NUTS

def model(tensor):
    # priors
    mu_b0 = pyro.sample('mu_b0', dist.Normal(0., 1.)).to(device)
    sig_b0 = pyro.sample('sig_b0', dist.HalfNormal(1.)).to(device)
    mu_b1 = pyro.sample('mu_b1', dist.Normal(0., 1.)).to(device)
    sig_b1 = pyro.sample('sig_b1', dist.HalfNormal(1.)).to(device)
    mu_b2 = pyro.sample('mu_b2', dist.Normal(0., 1.)).to(device)
    sig_b2 = pyro.sample('sig_b2', dist.HalfNormal(1.)).to(device)
    
    with pyro.plate('group', num_groups):
        b0 = pyro.sample('b0', dist.Normal(mu_b0, sig_b0)).to(device)
        b1 = pyro.sample('b1', dist.Normal(mu_b1, sig_b1)).to(device)
        b2 = pyro.sample('b2', dist.Normal(mu_b2, sig_b2)).to(device)
    
    mu = b0[tensor[:,0].long()] + b1[tensor[:,0].long()] * tensor[:,1] + b2[tensor[:,0].long()] * tensor[:,4]
    
    with pyro.plate('obs', len(tensor)):
        pyro.sample('hill_div_predict', dist.Normal(mu, 1.), obs=tensor[:,4])
        pyro.sample('species_count_predict', dist.Normal(mu, 1.), obs=tensor[:,5])

# Use GPU if available
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
pyro.set_rng_seed(1)
site = "OSBS"
df = gpd.read_file(f'data/{site}_gap_features.gpkg', driver='GPKG')

# remove geometry and turn geopandas into pandas
df = df.drop(columns=['geometry'])

df['group_index'] = df['sci_name'].astype('category').cat.codes
num_groups = len(df['group_index'].unique())

nuts_kernel = NUTS(model)
mcmc = MCMC(nuts_kernel, num_samples=1000, warmup_steps=200)

# Prepare tensor
tensor = torch.tensor(df[['group_index','dead_count','surface_dead','surface_alive', 'hill_div', 'species_count']].values).to(device)

mcmc.run(tensor)

# posterior
posterior = mcmc.get_samples

# inspect the posterior
print(posterior['mu_b0'].mean().item(), posterior['mu_b0'].std().item())
print(posterior['mu_b1'].mean().item(), posterior['mu_b1'].std().item())
print(posterior['mu_b2'].mean().item(), posterior['mu_b2'].std().item())

