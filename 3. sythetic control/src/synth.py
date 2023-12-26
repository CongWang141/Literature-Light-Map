import numpy as np
from sklearn.utils.validation import check_X_y, check_is_fitted, check_array
from sklearn.base import BaseEstimator, RegressorMixin

from scipy.optimize import minimize

class Synth(BaseEstimator, RegressorMixin):

    def __init__(self, ):
        pass
    
    def fit(self, data, outcome='packspercapita', id='state', time='year', treated='treated'):
        df = data.copy()
        df['treat'] = df.groupby(id)[treated].transform('max')
        df['post'] = df.groupby(time)[treated].transform('max')

        X = df[(df.post==0) & (df.treat==0)].pivot(index=time, columns=id, values=outcome)
        y = df[(df.post==0) & (df.treat==1)][outcome]
        X, y = check_X_y(X, y)

        initial_w = np.zeros(X.shape[1])

        def fun_obj(w, X, y):
            return np.sum((X @ w - y) ** 2)
        constraints = [{'type': 'eq', 'fun': lambda w: np.sum(w) - 1}]
        bounds = [(0, 1)] * X.shape[1]

        result = minimize(fun_obj, x0=initial_w, args=(X, y), method='SLSQP', bounds=bounds, constraints=constraints)

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