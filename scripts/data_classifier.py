"""Script to classify data using random forest"""

import paths
from classifier.random_forest import run_classifier

def main():
   scores, model = run_classifier(paths.training_data_path, paths.labelled_data)
   return None

if __name__ == '__main__':
    main()
