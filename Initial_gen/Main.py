from Generator import Generator
from XMLReader import XMLReader


def main():
    # Creating XMLReader object
    reader = XMLReader()
    reader.load()

    # Creating Generator object
    generator = Generator(reader)
    generator._build_configuration()
    generator._output_configuration()

__name__ == "__main__" and main()