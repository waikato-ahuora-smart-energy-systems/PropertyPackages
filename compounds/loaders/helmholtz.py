from loaders import register_loader

# Loaders get called when CompoundDB module is imported
@register_loader("chemsep")
def load():
    # Data is present within property package
    return "ye"
    