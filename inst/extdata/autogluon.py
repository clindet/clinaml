import sys
from autogluon.tabular import TabularDataset, TabularPredictor

m=sys.argv[1]
t=sys.argv[2]
o=sys.argv[3]

predictor = TabularPredictor.load(m)
testdat = TabularDataset(t)
predictor.predict_proba(testdat).to_csv(o)
