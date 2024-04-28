from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor

import numpy as np
from sklearn.base import BaseEstimator
from sklearn.tree import DecisionTreeClassifier


class RandomForestClassifierCustom(BaseEstimator):
    """
    Class to classify labels with Random Forest
    """
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

    def fit_tree(self, args) -> tuple:
        """
        Function to fit each tree of
        Random Forest
        Args:
            args: tuple of features, labels and
            number of estimators

        Returns:
            tuple of the tree and number
            of features for this tree

        """
        X, y, n_estim = args
        np.random.seed(seed=self.random_state + n_estim)
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

    def fit(self, X, y, n_jobs) -> BaseEstimator:
        """
        The function to fit Random Forest
        Args:
            X: features
            y: labels of interest
            n_jobs: the number of processes to run

        Returns:
            the fitted RandomForestClassifierCustom
            object
        """
        self.classes_ = sorted(np.unique(y))
        list_X = [X] * self.n_estimators
        list_y = [y] * self.n_estimators
        arguments = zip(list_X, list_y, range(self.n_estimators))
        with ProcessPoolExecutor(n_jobs) as pool:
            results = pool.map(self.fit_tree, arguments)
        for element in results:
            self.trees.append(element[0])
            self.feat_ids_by_tree.append(element[1])
        return self

    def predict_proba_tree(self, args) -> np.array:
        """
        Predict probabilities for each tree
        estimator in Random Forest
        Args:
            args: tuple of features, fitted trees and
            number of features for each tree

        Returns:
            array of predicted probabilities
             for each tree
        """
        X, tree, num_features = args
        prediction = tree.predict_proba(X[:, num_features])
        return prediction

    def predict_proba(self,
                      X,
                      n_jobs,
                      is_process=True) -> np.array:
        """
        Predict probabilities for the whole
        Random Forest
        Args:
            X: features
            n_jobs: the number of processes
            or flows
            is_process: parallel processes or flows
            By default is_process=True
            In case of is_process=False the run
            will be parallelized in flows

        Returns:
            array of predicted probabilities
            for Random Forest
        """
        list_X = [X] * self.n_estimators
        arguments = zip(list_X, self.trees, self.feat_ids_by_tree)
        if is_process:
            parallel_func = ProcessPoolExecutor
        else:
            parallel_func = ThreadPoolExecutor
        with parallel_func(n_jobs) as pool:
            predict_proba_list = pool.map(self.predict_proba_tree, arguments)
        return np.array(list(predict_proba_list)).mean(axis=0)

    def predict(self, X, n_jobs, is_process=True) -> np.ndarray:
        """
        According predicted probabilities predict
        labels
        Args:
            X: features
            n_jobs: the number of processes
            or flows
            is_process: parallel processes or flows
            By default is_process=True
            In case of is_process=False the run
            will be parallelized in flows

        Returns:
            array of obtained labels in
            classification

        """
        probas = self.predict_proba(X, n_jobs, is_process)
        predictions = np.argmax(probas, axis=1)
        return predictions
