import sklearn

def train(data, train_pct):
    X, y = train_test_split(X, y, train_pct)
    model = ...
    
    model.train(X,y)
    
    return model