
loaders_list = []

def loader(name):
    
    def decorator(func):
        loaders_list.append(func)
        return func
    
    return decorator