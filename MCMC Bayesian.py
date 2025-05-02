# ---------------- MCMC with additional Markov chain convergence diagnostics

import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt
from scipy.stats import norm
import seaborn as sns
from scipy.stats import beta, uniform
from scipy.stats import gaussian_kde
import xlwt
import arviz as az

# -------- Read data from Excel file
from openpyxl import load_workbook
# Open Excel file
workbook = load_workbook('random_numbers25.xlsx')  # Replace with your Excel file name

# Select the worksheet
sheet = workbook['1-1']

# Read data
data = []
for row in sheet.iter_rows(values_only=True):
    data.append(row)
# Close Excel file
workbook.close()


# Define the prior joint probability distribution (uniform distributions for a and b over [1,10])
def log_prior(a, b):
    return uniform.logpdf(a, 1, 9) + uniform.logpdf(b, 1, 9)  # uniform(lower bound, interval)


# Define the likelihood function
def log_likelihood(a, b, data):
    # Adjust normalization to avoid boundary issues
    data = np.array(data)
    eps = 0.000001
    normalized_data = (data - (0.3 - eps)) / (0.73 - 0.3 + 2 * eps)
    return np.sum(beta.logpdf(normalized_data, a, b))


# Define the posterior probability distribution
def log_posterior(a, b, data):
    return log_likelihood(a, b, data) + log_prior(a, b)

# Initialize parameters
a_current = 0.5 + np.random.randn()
b_current = 0.5 + np.random.rand()

# Store parameter values
a_values = []
b_values = []

# MCMC sampling for first chain
for i in range(10000):
    # Propose new parameter values from a normal distribution
    a_proposal = np.abs(np.random.normal(a_current, 0.4))
    b_proposal = np.abs(np.random.normal(b_current, 0.4))  # Ensure a and b remain positive

    # Calculate acceptance rate
    log_p_accept = min(0, log_posterior(a_proposal, b_proposal, data) - log_posterior(a_current, b_current, data))

    # Decide whether to accept new parameter values based on acceptance rate
    if np.log(np.random.rand()) < log_p_accept:
        a_current = a_proposal
        b_current = b_proposal
    a_values.append(a_current)
    b_values.append(b_current)

# Initialize second chain
a_current1 = 0.3 + np.random.randn()
b_current1 = 0.3 + np.random.rand()

# Store parameter values for second chain
a_values1 = []
b_values1 = []

# MCMC sampling for second chain
for i in range(10000):
    a_proposal1 = np.abs(np.random.normal(a_current1, 0.4))
    b_proposal1 = np.abs(np.random.normal(b_current1, 0.4))  # Ensure positivity

    log_p_accept = min(0, log_posterior(a_proposal1, b_proposal1, data) - log_posterior(a_current1, b_current1, data))

    if np.log(np.random.rand()) < log_p_accept:
        a_current1 = a_proposal1
        b_current1 = b_proposal1
    a_values1.append(a_current1)
    b_values1.append(b_current1)

# Combine values from both chains
a_values2 = np.concatenate([a_values, a_values1])
b_values2 = np.concatenate([b_values, b_values1])

print("Estimation of a:", np.mean(a_values2))
print("Estimation of b:", np.mean(b_values2))

# Trace plots for parameter 'a'
plt.figure(figsize=(8, 6))
plt.subplot(2, 1, 1)
plt.plot(a_values, color='lightblue', label='Chain 1 for a')
plt.plot(a_values1, color='salmon', label='Chain 2 for a')
plt.xlabel('Iteration')
plt.ylabel('a')
plt.title('Trace plot of a')
plt.legend()

# Trace plots for parameter 'b'
plt.subplot(2, 1, 2)
plt.plot(b_values, color='lightblue', label='Chain 1 for b')
plt.plot(b_values1, color='salmon', label='Chain 2 for b')
plt.xlabel('Iteration')
plt.ylabel('b')
plt.title('Trace plot of b')
plt.subplots_adjust(hspace=0.5)
plt.legend()
plt.show()

# Convert data to ArviZ format
data_dict = {
    "a": np.array([a_values, a_values1]),
    "b": np.array([b_values, b_values1])
}
idata = az.convert_to_inference_data(data_dict)

# Calculate Gelman-Rubin diagnostic
rhat_results = az.rhat(idata)

print("Gelman-Rubin Diagnostic (R-hat) for parameter a:", rhat_results['a'])
print("Gelman-Rubin Diagnostic (R-hat) for parameter b:", rhat_results['b'])

# Plot posterior distributions for a and b
plt.figure(figsize=(10, 12))

# Histogram and probability distribution for 'a'
plt.subplot(2, 1, 1)
n, bins, patches = plt.hist(a_values2, bins=50, alpha=0.6, color='tab:blue', edgecolor='black')
ax1 = plt.gca()
ax2 = ax1.twinx()
x1 = np.linspace(-0.5, 10.5, 1000)
pdf_uniform = uniform.pdf(x1, 0, 10)
ax2.plot(x1, pdf_uniform, color='tab:red')
ax2.set_ylabel('Probability')
ax1.set_ylabel('Frequency')

# Histogram and probability distribution for 'b'
plt.subplot(2, 1, 2)
n, bins, patches = plt.hist(b_values2, bins=30, alpha=0.6, color='tab:blue', edgecolor='black')
ax3 = plt.gca()
ax4 = ax3.twinx()
ax4.plot(x1, pdf_uniform, color='tab:red')
ax4.set_ylabel('Probability')
ax3.set_xlabel('b')
ax3.set_ylabel('Frequency')
plt.show()

# Plot prior and posterior distributions of variable X
x2_original = np.linspace(0.3, 0.73, 1000)
x2_normalized = (x2_original - 0.3) / 0.43
pdf_prior = beta.pdf(x2_normalized, 5, 5)
plt.plot(x2_original, pdf_prior, label='Prior distribution', color='tab:red')
pdf_posterior = beta.pdf(x2_normalized, np.mean(a_values2), np.mean(b_values2))
plt.plot(x2_original, pdf_posterior, label='Posterior distribution', color='tab:blue')
plt.xlabel('X')
plt.ylabel('Probability Density')
plt.grid(True, linestyle='--', linewidth=0.5)
plt.legend()
plt.show()
