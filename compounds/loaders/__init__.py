
loaders = {}

def loader(name):
    
    def decorator(func):
        loaders[name] = func
        return func
    
    return decorator

__all__ = ["loaders", "loader"]