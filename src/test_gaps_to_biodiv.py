
import torch
import pyro
import pyro.distributions as dist
import geopandas as gpd
import pandas as pd
from pyro.infer import MCMC, NUTS
import pyro
import pyro.distributions as dist
from pyro.infer import SVI, Trace_ELBO
from pyro.optim import Adam
import numpy as np
def model(tensor_data):
    # hyperparameters for normal priors
    mu_0, sigma_0 = 0., 1.

    # fixed effects coefficients
    beta_gap = pyro.sample("beta_gap", dist.Normal(mu_0, sigma_0))
    beta_density = pyro.sample("beta_density", dist.Normal(mu_0, sigma_0))
    beta_length = pyro.sample("beta_length", dist.Normal(mu_0, sigma_0))
    beta_ratio = pyro.sample("beta_ratio", dist.Normal(mu_0, sigma_0))
    
    # random effect coefficients
    # Define random effects for group1
    with pyro.plate("species_plate"):
        gamma_species = pyro.sample("gamma_species", dist.Normal(0, 1))
    
    # Define random effects for group2
    with pyro.plate("nlcd_plate"):
        gamma_nlcd = pyro.sample("gamma_nlcd", dist.Normal(0, 1))
    
    # Calculate linear combination for group1
    linear_combination_group1 = (tensor_data[:,2] + tensor_data[:,3] + tensor_data[:,4] + tensor_data[:,5]) * gamma_species
    
    # Calculate linear combination for group2
    linear_combination_group2 = (tensor_data[:,2] + tensor_data[:,3] + tensor_data[:,4] + tensor_data[:,5]) * gamma_nlcd
    
    # Model the relationship for both groups
    mu_group1 = beta_gap * tensor_data[:,2] + beta_density * tensor_data[:,3] + beta_length * tensor_data[:,4] + beta_ratio *tensor_data[:,5] + linear_combination_group1
    mu_group2 = beta_gap * tensor_data[:,2] + beta_density * tensor_data[:,3] + beta_length * tensor_data[:,4] + beta_ratio *tensor_data[:,5] + linear_combination_group2
   

    # Likelihood for both groups
    with pyro.plate("group1_data", len(tensor_data[:,0])):
        pyro.sample("obs_group1", dist.Normal(mu_group1, 0.1), obs=tensor_data[:,6][tensor_data[:,0].long()])
    
    with pyro.plate("group2_data", len(tensor_data[:,1])):
        pyro.sample("obs_group2", dist.Normal(mu_group2, 0.1), obs=tensor_data[:,6][tensor_data[:,1].long()])

# Define the guide
def guide(x):
    # Parameterized distribution for latent variables
    mu = pyro.param("mu", torch.tensor(0.0))
    sigma = pyro.param("sigma", torch.tensor(1.0), constraint=dist.constraints.positive)
    pyro.sample("latent_variable", dist.Normal(mu, sigma))

# Use GPU if available
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
pyro.set_rng_seed(1987)
site = "SERC"
# remove geometry and turn geopandas into pandas
df = gpd.read_file(f'data/{site}_gap_features.gpkg', driver='GPKG')
df = df.drop(columns=['geometry'])
df.to_csv(f'data/tmp{site}.csv')

df['focal_taxa'] = df['sci_name'].astype('category').cat.codes
df['ncdl'] = df['raster_value'].astype('category').cat.codes
df['dead_to_alive'] = df['surface_dead'] / df['surface_alive']
df['dead_to_alive'] = df['dead_to_alive'].fillna(0)
#fill inf with 1
df['dead_to_alive'] = df['dead_to_alive'].replace([np.inf, -np.inf], 1)
num_ecos = len(df['ncdl'].unique())
num_focal = len(df['focal_taxa'].unique())

# Prepare tensors
tensor_data = torch.tensor(df[['focal_taxa', 'ncdl', 'dead_count', 'network_length', 'surface_alive', 'surface_dead', 'hill_div']].values, dtype=torch.float).to(device)
#standardize each column of the tensor
for i in range(tensor_data.shape[1]):
    tensor_data[:,i] = (tensor_data[:,i] - tensor_data[:,i].mean()) / tensor_data[:,i].std()


#plot the distribution of each column of tensor
import matplotlib.pyplot as plt
#convert tenssurface
# write csv file of tmp



# Perform SVI
pyro.clear_param_store()
svi = SVI(model, guide, Adam({"lr": 0.01}), loss=Trace_ELBO())
num_iterations = 1000
for i in range(num_iterations):
    elbo = svi.step(tensor_data)

# After SVI training, obtain posterior samples using SVI's guide
posterior_samples = {name: pyro.param(name).detach().numpy() for name in pyro.get_param_store().get_all_param_names()}

# Assess convergence by plotting ELBO values

# Plot posterior distributions for each parameter
import matplotlib.pyplot as plt
import seaborn as sns

for param_name, samples in posterior_samples.items():
    plt.figure()
    sns.histplot(samples, kde=True)
    plt.title(f"Posterior Distribution of {param_name}")

# Compute credible intervals

# Calculate the median of each parameter's posterior distribution
medians = {param_name: np.median(samples) for param_name, samples in posterior_samples.items()}

# Perform hypothesis tests if needed

# Use convergence diagnostics if MCMC is used

# Interpret results and make inferences based on the analyses


