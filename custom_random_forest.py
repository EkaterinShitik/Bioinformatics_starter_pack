from multiprocessing import Pool

import numpy as np
from sklearn.base import BaseEstimator
from sklearn.tree import DecisionTreeClassifier


class RandomForestClassifierCustom(BaseEstimator):
    def __init__(self, n_estimators=10,
                 max_depth=None,
                 max_features=None,
                 random_state=None):
        self.n_estimators = n_estimators
        self.max_depth = max_depth
        self.max_features = max_features
        self.random_state = random_state
        self.trees = []
        self.feat_ids_by_tree = []
        self.classes_ = None

    def fit_tree(self, args):
        X, y, i = args
        np.random.seed(seed=self.random_state + i)
        rand_features = np.random.choice(X.shape[1],
                                         self.max_features,
                                         replace=False)
        rand_samples = np.random.choice(X.shape[0],
                                        X.shape[0],
                                        replace=True)
        X_sample = X[rand_samples][:, rand_features]
        y_sample = y[rand_samples]
        clf = DecisionTreeClassifier(max_depth=self.max_depth,
                                     max_features=self.max_features,
                                     random_state=self.random_state)
        clf.fit(X_sample, y_sample)
        return clf, rand_features

    def fit(self, X, y, n_jobs):
        self.classes_ = sorted(np.unique(y))
        list_X = [X] * self.n_estimators
        list_y = [y] * self.n_estimators
        arguments = zip(list_X, list_y, range(self.n_estimators))
        with Pool(n_jobs) as pool:
            results = pool.map(self.fit_tree, arguments)
        for element in results:
            self.trees.append(element[0])
            self.feat_ids_by_tree.append(element[1])
        return self

    def predict_proba(self, X, n_jobs):
        predict_proba_list = []
        number_tree = 0
        for tree in self.trees:
            prediction = tree.predict_proba(X[:, self.feat_ids_by_tree[number_tree]])
            predict_proba_list += [prediction]
            number_tree += 1
        return np.array(predict_proba_list).mean(axis=0)

    def predict(self, X, n_jobs):
        probas = self.predict_proba(X, n_jobs)
        predictions = np.argmax(probas, axis=1)

        return predictions
