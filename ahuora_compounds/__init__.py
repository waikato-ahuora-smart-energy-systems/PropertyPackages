import warnings

def deprecated(name):

    def decorator(func):

        def wrapper(*args, **kwargs):
            warnings.warn(f"{name} is deprecated and will be removed in future versions.", DeprecationWarning, stacklevel=2)
            return func(*args, **kwargs)
        
        return wrapper
    
    return decorator

__all__ = ["deprecated"]
