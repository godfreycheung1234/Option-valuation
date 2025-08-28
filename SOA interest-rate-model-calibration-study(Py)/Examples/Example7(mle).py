import pandas as pd
import numpy as np
from sklearn.linear_model import LinearRegression

file_path = r"\SOA paper\interest-calibration-1-120\InterestCalibrationv1\data\VeronesiTable14p7q5.csv"

data = pd.read_csv(file_path)
rt = data["rt"]
Delta = 1/252
N = len(rt)
y = np.array(rt[1:N])
x = np.array(rt[0:(N-1)])

X = x.reshape(-1, 1)
model = LinearRegression().fit(X, y)

y_pred = model.predict(X)
residuals = y - y_pred
sigma = np.sqrt(np.sum(residuals ** 2) / (N - 2))


mle_gamma = -np.log(model.coef_[0])/Delta
mle_rbar = model.intercept_/(1-model.coef_[0])
mle_sigma = sigma* ((2*mle_gamma/(1-model.coef_[0]**2))**0.5)

print(f"mle_gamma: {mle_gamma:.5f}")
print(f"mle_rbar: {mle_rbar:.5f}")
print(f"mle_sigma: {mle_sigma:.5f}")
