
import torch
import pyro
import pyro.distributions as dist
import geopandas as gpd
import pandas as pd
from pyro.infer import MCMC, NUTS


def model(tensor_data):
    # hyperparameters for normal priors
    mu_0, sigma_0 = 0., 1.
    tensor_num_focal = torch.tensor([num_focal]).to(device)
    tensor_num_ecos = torch.tensor([num_ecos]).to(device) 
    # fixed effects coefficients
    beta_gap = pyro.sample("beta_gap", dist.Normal(mu_0, sigma_0))
    beta_density = pyro.sample("beta_density", dist.Normal(mu_0, sigma_0))
    beta_length = pyro.sample("beta_length", dist.Normal(mu_0, sigma_0))
    beta_ratio = pyro.sample("beta_ratio", dist.Normal(mu_0, sigma_0))
    
    # random effect coefficients
    mu_species = pyro.sample("mu_species", dist.Normal(mu_0, sigma_0))
    sigma_species = pyro.sample("sigma_species", dist.HalfCauchy(1.))
    eta_species = pyro.sample("eta_species", dist.Normal(0., 1.), infer={"num_samples": tensor_num_focal})
    gamma_species = mu_species + sigma_species * eta_species # random intercept for species
    
    # random effect coefficients for nncdf
    mu_ncdf = pyro.sample("mu_ncdf", dist.Normal(mu_0, sigma_0))
    sigma_ncdf = pyro.sample("sigma_ncdf", dist.HalfCauchy(1.))
    eta_ncdf = pyro.sample("eta_ncdf", dist.Normal(0., 1.), infer={"num_samples": tensor_num_ecos})
    gamma_ncdf = mu_ncdf + sigma_ncdf * eta_ncdf  # random intercept for nncdf

    # linear predictor
    #    mu = beta_gap * size_gap + beta_length * network_length + beta_ratio * dead_alive_ratio + gamma_species[species]
    #mu = beta_gap * data[:,1] + beta_length * data[:,2] + beta_ratio * data[:,3] + gamma_species[data[:,0].long()]
    mu = (beta_gap * tensor_data[:,2] + beta_length * tensor_data[:,3] + beta_ratio * tensor_data[:,4] + beta_density * tensor_data[:,5]) * (gamma_species[tensor_data[:,0].long()] + gamma_ncdf[tensor_data[:,1].long()])

    # likelihood
    with pyro.plate("data", len(tensor_data)):
        pyro.sample("obs", dist.Normal(mu, 1.), obs=tensor_data[:,6])


# Use GPU if available
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
pyro.set_rng_seed(1987)
site = "HARV"
df = gpd.read_file(f'data/{site}_gap_features.gpkg', driver='GPKG')

# remove geometry and turn geopandas into pandas
df = df.drop(columns=['geometry'])
df['focal_taxa'] = df['sci_name'].astype('category').cat.codes
df['ncdl'] = df['raster_value'].astype('category').cat.codes
df['dead_to_alive'] = df['surface_dead'] / df['surface_alive']
num_ecos = len(df['ncdl'].unique())
num_focal = len(df['focal_taxa'].unique())

nuts_kernel = NUTS(model)
mcmc = MCMC(nuts_kernel, num_samples=1000, warmup_steps=200)

# Prepare tensors
tensor_data = torch.tensor(df[['focal_taxa', 'ncdl', 'dead_count', 'network_length', 'dead_to_alive', 'surface_dead', 'hill_div']].values, dtype=torch.float).to(device)

mcmc.run(tensor_data)

# posterior
posterior = mcmc.get_samples()
# save  mcmc  and posterior
torch.save(mcmc, f'models/{site}_mcmc.pt')
torch.save(nuts_kernel, f'models/{site}_model.pt')
torch.save(posterior, f'models/{site}_posterior.pt')

# load posterior
nuts_kernel = torch.load(f'models/{site}_model.pt')
mcmc = torch.load(f'models/{site}_mcmc.pt')
posterior = torch.load(f'models/{site}_posterior.pt')

# inspect the posterior
print(posterior['beta_gap'].mean().item(), posterior['beta_gap'].std().item())
print(posterior['beta_length'].mean().item(), posterior['beta_length'].std().item())
print(posterior['beta_ratio'].mean().item(), posterior['beta_ratio'].std().item())

