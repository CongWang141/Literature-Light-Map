import numpy as np
from sklearn.utils.validation import check_X_y, check_is_fitted, check_array
from sklearn.base import BaseEstimator, RegressorMixin

from scipy.optimize import minimize, Bounds, LinearConstraint

class Synth(BaseEstimator, RegressorMixin):

    def __init__(self, ):
        pass
    
    def fit(self, data, outcome='packspercapita', id='state', time='year', treated='treated'):
        df = data.copy()
        df['treat'] = df.groupby(id)[treated].transform('max')
        df['post'] = df.groupby(time)[treated].transform('max')

        X = df[(df.post==0) & (df.treat==0)].pivot(index=time, columns=id, values=outcome)
        y = df[(df.post==0) & (df.treat==1)].groupby(time)[outcome].mean()
        X, y = check_X_y(X, y)

        initial_w = np.ones(X.shape[1])/X.shape[1]
        v = np.diag(np.ones(X.shape[0])/X.shape[0])

        def fun_obj(w, X, y, v):
            return np.mean(np.sqrt((y - X @ w).T @ v @ (y - X @ w)))
        
        bounds = Bounds(lb=0, ub=1)
        constraints = LinearConstraint(np.ones(X.shape[1]), lb= 1, ub= 1)

        result = minimize(fun_obj, x0=initial_w, args=(X, y, v), method='SLSQP', bounds=bounds, constraints=constraints)

        self.X_ = X
        self.y_ = y
        self.w_ = result.x
        self.df_ = df

        self.id_ = id
        self.outcome_ = outcome
        self.time_ = time
        self.treated_ = treated

        self.is_fitted_ = True
        return self
    
    def predict(self):
        check_is_fitted(self)
        X = self.df_[self.df_.treat==0].pivot(index=self.time_, columns=self.id_, values=self.outcome_)
        X = check_array(X)
        return X @ self.w_