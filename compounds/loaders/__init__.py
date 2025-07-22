
loaders_list = []

print("Compounds loaders initialized.")

def loader(name):

    print("Loader registered:", name)
    
    def decorator(func):
        loaders_list.append(func)
        return func
    
    return decorator