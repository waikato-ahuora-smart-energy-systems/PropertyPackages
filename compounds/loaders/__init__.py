loaders = {}

def register_loader(name):
    """
    Decorator to register a loader function for a specific property package.
    
    Args:
        name (str): The name of the property package to register the loader for.
    
    Returns:
        function: The decorated loader function.
    """
    def decorator(func):
        loaders[name] = func
        return func
    
    return decorator